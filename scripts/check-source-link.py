#!/usr/bin/env python3
"""Verify committed source retrieval against a matching portable-PDB manifest.

Embedded documents are disclosed, but this gate does not extract or verify their
content. The metadata tool verifies the DLL/PDB identity before writing a manifest.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from collections.abc import Callable, Mapping
from pathlib import Path
from urllib.error import URLError
from urllib.parse import quote, unquote, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

MAX_MANIFEST_BYTES = 4 * 1024 * 1024
MAX_SOURCE_BYTES = 8 * 1024 * 1024
MAX_TOTAL_SOURCE_BYTES = 64 * 1024 * 1024
MAX_DOCUMENTS = 1000
IO_TIMEOUT_SECONDS = 10
DOWNLOAD_DEADLINE_SECONDS = 30
HASH_ALGORITHMS = {
    "ff1816ec-aa5e-4d10-87f7-6f4963833460": ("sha1", 40),
    "8829d00f-11b8-4213-878b-770e8597ac16": ("sha256", 64),
}


class RefuseRedirects(HTTPRedirectHandler):
    """Keep source retrieval on the validated host and exact committed URL."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError("Source retrieval redirects are prohibited.")


def download(url: str) -> bytes:
    """Read bounded bytes with verified HTTPS, idle timeout, and an elapsed guard.

    read1 performs one buffered read, allowing the elapsed guard to run between
    socket operations. The elapsed limit may be exceeded by one idle timeout;
    DNS lookup retains the host resolver's own timeout policy.
    """
    deadline = time.monotonic() + DOWNLOAD_DEADLINE_SECONDS
    request = Request(url, headers={"User-Agent": "DotNetMatrix-SourceLink-Verification"})
    with build_opener(RefuseRedirects()).open(request, timeout=IO_TIMEOUT_SECONDS) as response:
        if response.status != 200 or response.geturl() != url:
            raise ValueError("Source retrieval did not return the exact requested URL with HTTP 200.")
        content_length = response.headers.get("Content-Length")
        if content_length is not None and not 0 <= int(content_length) <= MAX_SOURCE_BYTES:
            raise ValueError("Source response exceeds its byte limit.")
        chunks: list[bytes] = []
        length = 0
        while True:
            if time.monotonic() >= deadline:
                raise TimeoutError("Source response exceeded its elapsed-time limit.")
            chunk = response.read1(min(64 * 1024, MAX_SOURCE_BYTES - length + 1))
            if not chunk:
                return b"".join(chunks)
            length += len(chunk)
            if length > MAX_SOURCE_BYTES:
                raise ValueError("Source response exceeds its byte limit.")
            chunks.append(chunk)


def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def load_manifest(path: Path) -> object:
    with path.open("rb") as stream:
        payload = stream.read(MAX_MANIFEST_BYTES + 1)
    if len(payload) > MAX_MANIFEST_BYTES:
        raise ValueError("Source Link manifest exceeds its byte limit.")
    return json.loads(payload, object_pairs_hook=unique_object)


def validate_url_template(template: object, commit: str) -> str:
    if not isinstance(template, str):
        raise ValueError("Source Link mapping URL must be a string.")
    parsed = urlsplit(template)
    if parsed.scheme != "https" or parsed.netloc != "raw.githubusercontent.com" or parsed.query or parsed.fragment:
        raise ValueError("Source Link mapping must use the approved GitHub raw HTTPS host without query or fragment.")
    match = re.fullmatch(r"/firestrand/DotNetMatrix/([0-9a-fA-F]{40})/(.+)", parsed.path)
    if match is None:
        raise ValueError("Source Link mapping must identify the DotNetMatrix repository and a full committed SHA.")
    if match[1].lower() != commit:
        raise ValueError("Source Link SHA drift: mapping does not identify the expected commit.")
    relative = unquote(match[2])
    if any(part in {"", ".", ".."} for part in relative.split("/")) or "\\" in relative:
        raise ValueError("Source Link mapping contains an invalid repository path.")
    if any(ord(character) < 32 for character in relative) or relative.count("*") > 1:
        raise ValueError("Source Link mapping contains invalid characters or wildcards.")
    if "*" in relative and not relative.endswith("*"):
        raise ValueError("Source Link wildcard must end its URL template.")
    return template


def source_url(name: str, mappings: Mapping[str, str], commit: str) -> str:
    matches: list[str] = []
    for pattern, template in mappings.items():
        if pattern.endswith("*") and name.startswith(pattern[:-1]):
            suffix = name[len(pattern) - 1:].replace("\\", "/")
            if not suffix or any(part in {"", ".", ".."} for part in suffix.split("/")):
                raise ValueError("Document mapping substitution contains an invalid repository path.")
            matches.append(template[:-1] + quote(suffix, safe="/"))
        elif pattern == name:
            matches.append(template)
    if len(matches) != 1:
        raise ValueError(f"Document must have exactly one Source Link mapping: {name}")
    return validate_url_template(matches[0], commit)


def check_manifest(
    manifest: object,
    expected_commit: str,
    downloader: Callable[[str], bytes] = download,
) -> dict[str, object]:
    if re.fullmatch(r"[0-9a-fA-F]{40}", expected_commit) is None:
        raise ValueError("Expected commit must be a full 40-character hexadecimal Git SHA.")
    commit = expected_commit.lower()
    if not isinstance(manifest, dict) or set(manifest) != {"schemaVersion", "sourceLink", "documents"}:
        raise ValueError("Invalid Source Link manifest schema.")
    if type(manifest["schemaVersion"]) is not int or manifest["schemaVersion"] != 1:
        raise ValueError("Unknown Source Link manifest schema version.")
    source_link = manifest["sourceLink"]
    if not isinstance(source_link, dict) or set(source_link) != {"documents"}:
        raise ValueError("Missing or invalid Source Link map.")
    mappings = source_link["documents"]
    if not isinstance(mappings, dict) or not mappings:
        raise ValueError("Source Link mappings must not be empty.")
    for pattern, template in mappings.items():
        if not isinstance(pattern, str) or not pattern or any(ord(character) < 32 for character in pattern):
            raise ValueError("Source Link mapping path must be a nonempty string without control characters.")
        validate_url_template(template, commit)
        if pattern.count("*") > 1 or ("*" in pattern and not pattern.endswith("*")):
            raise ValueError("Source Link mapping wildcard must occur once at the end.")
        if pattern.endswith("*") != template.endswith("*"):
            raise ValueError("Source Link mapping path and URL wildcards must agree.")

    documents = manifest["documents"]
    if not isinstance(documents, list) or not 1 <= len(documents) <= MAX_DOCUMENTS:
        raise ValueError("Source Link documents must be nonempty and within the count limit.")
    names: set[str] = set()
    embedded: list[str] = []
    planned: list[tuple[str, str, str, str]] = []
    for document in documents:
        if not isinstance(document, dict) or set(document) != {"name", "hashAlgorithm", "hash", "embeddedSource"}:
            raise ValueError("Invalid Source Link document schema.")
        name = document["name"]
        if not isinstance(name, str) or not name or any(ord(character) < 32 for character in name):
            raise ValueError("Document name must be nonempty and contain no control characters.")
        if name in names:
            raise ValueError(f"Duplicate document: {name}")
        names.add(name)
        algorithm = document["hashAlgorithm"]
        if not isinstance(algorithm, str) or algorithm not in HASH_ALGORITHMS:
            raise ValueError(f"Unknown document hash algorithm: {name}")
        hash_name, hash_length = HASH_ALGORITHMS[algorithm]
        checksum = document["hash"]
        if not isinstance(checksum, str) or re.fullmatch(rf"[0-9a-fA-F]{{{hash_length}}}", checksum) is None:
            raise ValueError(f"Malformed document checksum: {name}")
        if type(document["embeddedSource"]) is not bool:
            raise ValueError(f"Embedded-source flag must be Boolean: {name}")
        if document["embeddedSource"]:
            embedded.append(name)
        else:
            planned.append((name, source_url(name, mappings, commit), hash_name, checksum.lower()))

    verified: list[str] = []
    total_bytes = 0
    for name, url, algorithm, expected_hash in planned:
        try:
            content = downloader(url)
        except (OSError, URLError, ValueError) as error:
            raise ValueError(f"Source download failed for {name}: {error}") from error
        if not isinstance(content, bytes) or len(content) > MAX_SOURCE_BYTES:
            raise ValueError(f"Source download has invalid type or exceeds its byte limit: {name}")
        total_bytes += len(content)
        if total_bytes > MAX_TOTAL_SOURCE_BYTES:
            raise ValueError("Source retrieval exceeds its total byte limit.")
        actual_hash = hashlib.new(algorithm, content).hexdigest()
        if actual_hash != expected_hash:
            raise ValueError(f"Source checksum mismatch: {name}")
        verified.append(name)
    return {
        "expectedCommit": commit,
        "verifiedSourceDocuments": verified,
        "embeddedDocumentsNotVerified": embedded,
        "embeddedSourcePolicy": "Disclosed; network retrieval skipped; embedded content was not verified by this gate.",
    }


def check(path: Path, expected_commit: str, downloader: Callable[[str], bytes] = download) -> dict[str, object]:
    return check_manifest(load_manifest(path), expected_commit, downloader)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("expected_commit")
    arguments = parser.parse_args()
    try:
        report = check(arguments.manifest, arguments.expected_commit)
    except (ValueError, OSError) as error:
        print(f"Source Link verification failed: {error}", file=sys.stderr)
        return 1
    sys.stdout.buffer.write((json.dumps(report, indent=2) + "\n").encode("utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
