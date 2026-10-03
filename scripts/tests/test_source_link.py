"""Test Source Link trust boundaries without downloading fixture source."""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "source_link_gate", Path(__file__).resolve().parents[1] / "check-source-link.py"
)
assert SPEC is not None and SPEC.loader is not None
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)

COMMIT = "a" * 40
BASE_URL = f"https://raw.githubusercontent.com/firestrand/DotNetMatrix/{COMMIT}/"
SOURCE = b"namespace DotNetMatrix;\n"
SHA256_GUID = "8829d00f-11b8-4213-878b-770e8597ac16"


class SourceLinkTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = {
            "schemaVersion": 1,
            "sourceLink": {"documents": {"/_/*": BASE_URL + "*"}},
            "documents": [{
                "name": "/_/DotNetMatrix/Matrix.cs",
                "hashAlgorithm": SHA256_GUID,
                "hash": hashlib.sha256(SOURCE).hexdigest().upper(),
                "embeddedSource": False,
            }],
        }
        self.urls: list[str] = []

    def fetch(self, url: str) -> bytes:
        self.urls.append(url)
        return SOURCE

    def verify(self) -> dict[str, object]:
        return GATE.check_manifest(self.manifest, COMMIT, self.fetch)

    def test_wildcard_substitution_checksum_and_commit_binding(self) -> None:
        report = self.verify()
        self.assertEqual([BASE_URL + "DotNetMatrix/Matrix.cs"], self.urls)
        self.assertEqual(COMMIT, report["expectedCommit"])
        self.assertEqual(["/_/DotNetMatrix/Matrix.cs"], report["verifiedSourceDocuments"])
        self.assertEqual([], report["embeddedDocumentsNotVerified"])

    def test_exact_mapping_and_percent_encoded_unicode_substitution(self) -> None:
        self.manifest["sourceLink"]["documents"] = {
            "/_/DotNetMatrix/Matrix.cs": BASE_URL + "DotNetMatrix/Matrix.cs"
        }
        self.verify()
        self.assertEqual([BASE_URL + "DotNetMatrix/Matrix.cs"], self.urls)
        self.manifest["sourceLink"]["documents"] = {"/_/*": BASE_URL + "*"}
        self.manifest["documents"][0]["name"] = "/_/My source/é.cs"
        self.verify()
        self.assertEqual(BASE_URL + "My%20source/%C3%A9.cs", self.urls[-1])

    def test_windows_document_path_substitution(self) -> None:
        self.manifest["sourceLink"]["documents"] = {"C:\\work\\*": BASE_URL + "*"}
        self.manifest["documents"][0]["name"] = "C:\\work\\DotNetMatrix\\Matrix.cs"
        self.verify()
        self.assertEqual([BASE_URL + "DotNetMatrix/Matrix.cs"], self.urls)

    def test_embedded_documents_are_disclosed_without_claiming_verification(self) -> None:
        embedded = copy.deepcopy(self.manifest["documents"][0])
        embedded.update(name="/generated/AssemblyInfo.cs", embeddedSource=True)
        self.manifest["documents"].append(embedded)
        report = self.verify()
        self.assertEqual(1, len(self.urls))
        self.assertEqual([embedded["name"]], report["embeddedDocumentsNotVerified"])
        self.assertIn("not verified", report["embeddedSourcePolicy"])
        self.assertNotIn(embedded["name"], report["verifiedSourceDocuments"])

    def test_unknown_or_unapproved_urls_rejected_before_download(self) -> None:
        invalid_urls = [
            BASE_URL.replace("https:", "http:"),
            BASE_URL.replace("raw.githubusercontent.com", "example.invalid"),
            BASE_URL.replace("raw.githubusercontent.com", "raw.githubusercontent.com.evil.invalid"),
            BASE_URL.replace("raw.githubusercontent.com", "user@raw.githubusercontent.com"),
            BASE_URL.replace("raw.githubusercontent.com", "raw.githubusercontent.com:443"),
            BASE_URL.replace("firestrand/DotNetMatrix", "other/Repository"),
            BASE_URL + "Matrix.cs?token=value",
            BASE_URL + "Matrix.cs#fragment",
            BASE_URL.replace(COMMIT, "main"),
            BASE_URL + "../Matrix.cs",
            BASE_URL + "%2e%2e/Matrix.cs",
            BASE_URL + "sub//Matrix.cs",
        ]
        for url in invalid_urls:
            with self.subTest(url=url):
                self.manifest["sourceLink"]["documents"] = {"/_/*": url + "*"}
                with self.assertRaises(ValueError):
                    self.verify()
        self.assertEqual([], self.urls)

    def test_expected_sha_drift_or_non_full_commit_rejected(self) -> None:
        for commit in ["b" * 40, COMMIT[:7], "main", "g" * 40, ""]:
            with self.subTest(commit=commit), self.assertRaisesRegex(ValueError, "commit|SHA"):
                GATE.check_manifest(self.manifest, commit, self.fetch)
        self.assertEqual([], self.urls)

    def test_hash_mismatch_and_unknown_or_malformed_hash_rejected(self) -> None:
        self.manifest["documents"][0]["hash"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.verify()
        for field, value in [("hash", "invalid"), ("hashAlgorithm", "unknown"), ("hashAlgorithm", None)]:
            with self.subTest(field=field):
                self.manifest["documents"][0][field] = value
                with self.assertRaises(ValueError):
                    self.verify()

    def test_known_sha1_portable_pdb_algorithm(self) -> None:
        self.manifest["documents"][0].update(
            hashAlgorithm="ff1816ec-aa5e-4d10-87f7-6f4963833460", hash=hashlib.sha1(SOURCE).hexdigest()
        )
        self.assertEqual(1, len(self.verify()["verifiedSourceDocuments"]))

    def test_duplicate_empty_unmapped_or_ambiguous_documents_rejected(self) -> None:
        original = copy.deepcopy(self.manifest)
        for documents in [[], [original["documents"][0]] * 2, [{**original["documents"][0], "name": ""}],
                          [{**original["documents"][0], "name": "/unknown/Matrix.cs"}]]:
            with self.subTest(documents=documents):
                self.manifest = {**original, "documents": documents}
                with self.assertRaises(ValueError):
                    self.verify()
        self.manifest = original
        self.manifest["sourceLink"]["documents"]["/_/DotNetMatrix/*"] = BASE_URL + "DotNetMatrix/*"
        with self.assertRaisesRegex(ValueError, "exactly one"):
            self.verify()
        self.assertEqual([], self.urls)

    def test_traversal_substitution_or_inconsistent_wildcards_rejected(self) -> None:
        self.manifest["documents"][0]["name"] = "/_/../private.cs"
        with self.assertRaisesRegex(ValueError, "invalid repository path"):
            self.verify()
        self.manifest["documents"][0]["name"] = "/_/DotNetMatrix/Matrix.cs"
        for pattern, url in [("/_/*", BASE_URL + "Matrix.cs"), ("/_/Matrix.cs", BASE_URL + "*"),
                             ("/_/*/*", BASE_URL + "*"), ("/_/*", BASE_URL + "*/Matrix.cs")]:
            self.manifest["sourceLink"]["documents"] = {pattern: url}
            with self.subTest(pattern=pattern, url=url), self.assertRaises(ValueError):
                self.verify()
        self.assertEqual([], self.urls)

    def test_invalid_schema_flags_and_maps_rejected_before_download(self) -> None:
        original = copy.deepcopy(self.manifest)
        malformed = [None, {}, {**original, "schemaVersion": True}, {**original, "schemaVersion": 2},
                     {**original, "sourceLink": None}, {**original, "sourceLink": {"documents": {}}},
                     {**original, "extra": True}, {**original, "documents": [{**original["documents"][0], "embeddedSource": 1}]}]
        for manifest in malformed:
            with self.subTest(manifest=manifest), self.assertRaises(ValueError):
                GATE.check_manifest(manifest, COMMIT, self.fetch)
        self.assertEqual([], self.urls)

    def test_download_failure_and_oversized_download_fail_closed(self) -> None:
        def unavailable(url: str) -> bytes:
            raise OSError("Synthetic network failure")

        with self.assertRaisesRegex(ValueError, "Source download failed"):
            GATE.check_manifest(self.manifest, COMMIT, unavailable)
        with patch.object(GATE, "MAX_SOURCE_BYTES", 3), self.assertRaisesRegex(ValueError, "byte limit"):
            self.verify()
        with patch.object(GATE, "MAX_TOTAL_SOURCE_BYTES", 3), self.assertRaisesRegex(ValueError, "total byte limit"):
            self.verify()

    def test_manifest_file_boundary_rejects_duplicate_keys_and_oversize(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps(self.manifest), encoding="utf-8")
            self.assertEqual(COMMIT, GATE.check(path, COMMIT, self.fetch)["expectedCommit"])
            path.write_text('{"schemaVersion":1,"schemaVersion":1}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate JSON key"):
                GATE.load_manifest(path)
            with patch.object(GATE, "MAX_MANIFEST_BYTES", 2), self.assertRaisesRegex(ValueError, "byte limit"):
                GATE.load_manifest(path)

    def test_http_adapter_enforces_size_exact_url_timeout_and_no_redirect(self) -> None:
        class Response(io.BytesIO):
            status = 200
            headers = {}

            def geturl(self) -> str:
                return BASE_URL + "Matrix.cs"

        url = BASE_URL + "Matrix.cs"
        with patch.object(GATE, "build_opener") as opener:
            opener.return_value.open.return_value = Response(SOURCE)
            self.assertEqual(SOURCE, GATE.download(url))
            self.assertEqual(GATE.IO_TIMEOUT_SECONDS, opener.return_value.open.call_args.kwargs["timeout"])
            response = Response(SOURCE)
            response.headers = {"Content-Length": str(GATE.MAX_SOURCE_BYTES + 1)}
            opener.return_value.open.return_value = response
            with self.assertRaisesRegex(ValueError, "byte limit"):
                GATE.download(url)
            opener.return_value.open.return_value = Response(SOURCE)
            with self.assertRaisesRegex(ValueError, "exact requested URL"):
                GATE.download(BASE_URL + "Other.cs")
            opener.return_value.open.return_value = Response(SOURCE)
            with patch.object(GATE.time, "monotonic", side_effect=[0, 31]), self.assertRaises(TimeoutError):
                GATE.download(url)
        with self.assertRaisesRegex(ValueError, "redirects"):
            GATE.RefuseRedirects().redirect_request(None, None, 302, "redirect", {}, "https://example.invalid")


if __name__ == "__main__":
    unittest.main()
