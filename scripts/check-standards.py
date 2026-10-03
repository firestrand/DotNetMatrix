#!/usr/bin/env python3
"""Offline standards lint; exception authorization comes only from protected evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

STANDARD_VERSION = "1.1.0-draft.1"
ADOPTED_SOURCE_SHA256 = "8c56930ec1cd26e8b58a273d47cb0b9713f55e6d4961b6d29a636ffaf4f0ebe2"
RULE_PATTERN = re.compile(r"DOTNET-[A-Z]+-[0-9]{3}")
MARKER = re.compile(r"OVERRIDE\((DOTNET-[A-Z]+-[0-9]{3}),\s*(EXC-[0-9]{4}-[0-9]{3})\)")
EXCEPTION_PATTERN = re.compile(r"EXC-[0-9]{4}-[0-9]{3}\Z")
WARNING_CLAUSE = "CI and release publication **MUST** fail on the configured compiler/analyzer warnings."


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_json(path: Path):
    # JSON is a YAML 1.2 subset. Only this documented subset is accepted for the registry.
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate object key {key!r}")
            result[key] = value
        return result
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique)


def source_rules(source: str) -> list[dict]:
    lines = source.splitlines()
    definitions = []
    for index, line in enumerate(lines):
        group = line.startswith("**Rule group ")
        individual = line.startswith("**Standard `DOTNET-") or line.startswith("**`DOTNET-")
        if group or individual:
            match = RULE_PATTERN.search(line)
            if match:
                end = index + 1
                if group:
                    while end < len(lines) and not lines[end].startswith(("**Rule group ", "## ")):
                        end += 1
                title = line.split(" — ", 1)[-1].strip("* .") if group else line.strip("*")
                definitions.append({"id": match.group(), "title": title, "kind": "group" if group else "rule",
                                    "source": {"path": "docs/standards/engineering-standard.md", "start_line": index + 1, "end_line": end}})
    return definitions


def exact_normative_sentence(clause: str, bounded: str) -> bool:
    if not re.search(r"\*\*(?:MUST|MUST NOT|SHOULD|SHOULD NOT)(?::)?\*\*", clause) or not clause.endswith((".", "!", "?")):
        return False
    # Exact whole sentence, not a selectively quoted fragment of the requirement.
    normalized = re.sub(r"(?m)^\*\*[^*\n]+[.:]\*\*\s*", "", bounded)
    normalized = re.sub(r"[ \t]+", " ", normalized)
    for match in re.finditer(re.escape(clause), normalized):
        before = normalized[:match.start()].rstrip(" \t")
        after = normalized[match.end():]
        if (not before or before.endswith((".", "!", "?", "\n"))) and (not after or after[0].isspace()):
            return True
    return False


def validate_catalog(root: Path) -> tuple[dict, list[str]]:
    errors = []
    catalog = load_json(root / "docs/standards/dotnet-rules.json")
    source_bytes = (root / "docs/standards/engineering-standard.md").read_bytes()
    if catalog.get("schema_version") != 1 or catalog.get("standard_version") != STANDARD_VERSION:
        errors.append("DOTNET-GOV-001: unsupported catalog schema/standard version")
    if catalog.get("source_sha256") != sha256(source_bytes) or sha256(source_bytes) != ADOPTED_SOURCE_SHA256:
        errors.append("DOTNET-GOV-001: adopted source SHA-256 mismatch")
    expected = source_rules(source_bytes.decode("utf-8"))
    entries = catalog.get("rules", [])
    if len(entries) != 43 or len({entry.get("id") for entry in entries}) != len(entries):
        errors.append("DOTNET-GOV-001: expected 43 unique catalog definitions")
    if [{key: entry.get(key) for key in ("id", "title", "kind", "source")} for entry in entries] != expected:
        errors.append("DOTNET-GOV-001: catalog IDs/titles/source ranges differ from adopted revision")
    for entry in entries:
        source_range = expected[next((index for index, rule in enumerate(expected) if rule["id"] == entry.get("id")), 0)]["source"]
        bounded = "\n".join(source_bytes.decode("utf-8").splitlines()[source_range["start_line"] - 1:source_range["end_line"]])
        clauses = entry.get("clauses", [])
        if not isinstance(clauses, list) or not all(isinstance(clause, str) for clause in clauses) or len(clauses) != len(set(clauses)):
            errors.append(f"{entry.get('id')}: malformed/duplicate clause selectors")
        else:
            for clause in clauses:
                named = f"**Clause `{clause}`:"
                if named not in bounded and not exact_normative_sentence(clause, bounded):
                    errors.append(f"{entry.get('id')}: clause selector is not a declared key or exact normative source sentence")
        if not entry.get("owner") or entry.get("status") != "adopted-for-project-audit":
            errors.append(f"{entry.get('id')}: missing owner role/adoption status")
        if not entry.get("enforcement") or any(not re.match(r"(?:gate:|test:|review:).+", item) for item in entry["enforcement"]):
            errors.append(f"{entry.get('id')}: missing gate/test/review enforcement contract")
    return catalog, errors


def utc_instant(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None or parsed.utcoffset() != timedelta(0):
        raise ValueError("evidence timestamp must explicitly be UTC")
    return parsed


def validate_record(record: dict, catalog: dict, now: datetime, registry_digest: str,
                    evidence: dict | None = None) -> list[str]:
    errors = []
    ident = record.get("id", "<missing>")
    rule = record.get("rule", "<missing>")
    prefix = f"{ident} {rule} {record.get('scope', {})}: "
    required = {"id", "standard_version", "rule", "clause", "diagnostics", "scope", "reason", "owner", "approved_by",
                "security_impact", "compensating_controls", "created_on", "expires_on", "remediation_issue", "status"}
    if not required.issubset(record):
        return [prefix + "missing registry fields: " + ", ".join(sorted(required - record.keys()))]
    if not isinstance(ident, str) or not EXCEPTION_PATTERN.fullmatch(ident):
        errors.append("malformed exception ID")
    if record["standard_version"] != STANDARD_VERSION:
        errors.append("unknown standard version")
    known = next((item for item in catalog["rules"] if item["id"] == rule), None)
    if known is None:
        errors.append("unknown rule ID")
    elif record["clause"] not in known.get("clauses", []):
        errors.append("unknown/ambiguous clause selector")
    if not isinstance(record["scope"], dict) or not record["scope"].get("files") or not record["scope"].get("symbols"):
        errors.append("exact file and symbol scope required")
    else:
        if not isinstance(record["scope"]["files"], list) or not isinstance(record["scope"]["symbols"], list) or not all(isinstance(symbol, str) and symbol for symbol in record["scope"]["symbols"]):
            errors.append("scope files/symbols must be nonempty string lists")
        for file in record["scope"]["files"]:
            if not isinstance(file, str) or Path(file).is_absolute() or ".." in Path(file).parts or any(char in file for char in "*?["):
                errors.append("scope must use exact repository-relative files")
    for field in ("reason", "owner", "security_impact", "compensating_controls", "diagnostics"):
        if not record[field]:
            errors.append(f"missing {field}")
    try:
        created = date.fromisoformat(record["created_on"])
        expiry = date.fromisoformat(record["expires_on"])
        if expiry <= created:
            errors.append("expiration must follow creation")
        if now.date() >= expiry:
            errors.append("expired at 00:00 UTC on expires_on")
    except (TypeError, ValueError):
        errors.append("malformed calendar date")
    if record["status"] != "active":
        errors.append(f"{record['status']} record cannot authorize a finding")
    if not record["approved_by"] or record["approved_by"] == record["owner"]:
        errors.append("missing independent approval")
    if not record["remediation_issue"] or not str(record["remediation_issue"]).startswith("https://"):
        errors.append("missing HTTPS remediation/review issue")
    if evidence is None:
        errors.append("missing protected external approval/issue evidence")
    else:
        if evidence.get("registry_sha256") != registry_digest or evidence.get("owner_review_verified") is not True:
            errors.append("registry revision/protected owner review not verified")
        approved = evidence.get("exceptions", {}).get(ident, {})
        if approved.get("record_sha256") != sha256(json.dumps(record, sort_keys=True).encode()):
            errors.append("scope/record changed after independent review")
        if approved.get("approved_by") != record["approved_by"] or approved.get("approval_verified") is not True:
            errors.append("independent approval not verified")
        if approved.get("issue") != record["remediation_issue"] or approved.get("issue_state") != "open":
            errors.append("issue closed, unavailable, or unknown")
        try:
            age = now - utc_instant(approved["checked_at"])
            if age < timedelta(0) or age > timedelta(hours=24):
                errors.append("issue evidence stale/future (maximum age 24 hours)")
        except (KeyError, TypeError, ValueError):
            errors.append("issue evidence missing/malformed/API failure")
    return [prefix + error for error in errors]


def lexical_source(text: str) -> str:
    """Remove ordinary comments/string bodies, preserving positions and newlines."""
    pattern = re.compile(r'//[^\n]*|/\*.*?\*/|@"(?:""|[^"])*"|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', re.DOTALL)
    return pattern.sub(lambda match: "".join("\n" if char == "\n" else " " for char in match.group()), text)


def symbol_body(lines: list[str], index: int) -> tuple[str | None, int | None]:
    """Conservative method-body detection; unsupported forms fail closed."""
    text = "\n".join(lines)
    clean = lexical_source(text)
    position = sum(len(line) + 1 for line in lines[:index])
    before = clean[:position]
    namespace = re.search(r"\bnamespace\s+([A-Za-z_][\w.]*)", before)
    method_pattern = r"^\s*(?:public|private|protected|internal)\s+(?:(?:static|virtual|override|async|sealed|new)\s+)*[\w<>?\[\].]+\s+([A-Za-z_]\w*)\s*\([^;{}]*\)\s*\{"
    methods = list(re.finditer(method_pattern, before, re.MULTILINE))
    if not namespace or not methods:
        return None, None
    method = methods[-1]
    classes = re.findall(r"\b(?:class|struct|record)\s+([A-Za-z_]\w*)", before[:method.start()])
    if not classes:
        return None, None
    # Raw strings are unsupported here rather than incorrectly expanding scope.
    if '\"\"\"' in text[method.start():position]:
        return None, None
    opening = method.end() - 1
    depth = 0
    closing = None
    for offset in range(opening, len(clean)):
        if clean[offset] == "{":
            depth += 1
        elif clean[offset] == "}":
            depth -= 1
            if depth == 0:
                closing = offset
                break
    if closing is None or not opening < position < closing:
        return None, None
    return namespace.group(1) + "." + classes[-1] + "." + method.group(1), clean[:closing].count("\n") + 1


def containing_symbol(lines: list[str], index: int) -> str | None:
    return symbol_body(lines, index)[0]


def source_paths(root: Path):
    excluded = {".git", ".nuget", ".serena", "obj", "bin", "artifacts", "BenchmarkDotNet.Artifacts"}
    for path in root.rglob("*"):
        if not path.is_file() or any(part in excluded for part in path.relative_to(root).parts):
            continue
        if path.suffix in {".cs", ".csproj", ".props", ".targets"} or path.name in {".editorconfig", ".globalconfig"}:
            yield path


def diagnostic_ids(text: str) -> list[str]:
    identifiers = re.findall(r"\b(?:(?:CS|CA|IDE|SYSLIB|NU)\d+|\d+)\b", text)
    return ["CS" + item if item.isdigit() else item for item in identifiers]


def findings(root: Path) -> tuple[list[dict], list[str]]:
    found, errors = [], []
    for path in source_paths(root):
        relative = path.relative_to(root).as_posix()
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()

        def add(index: int, diagnostics: list[str], body: bool = False):
            context = "\n".join(lines[max(0, index - 2):index + 1])
            marker = MARKER.search(context)
            symbol, body_end = symbol_body(lines, index) if body else ("configuration:" + relative, None)
            if body:
                remaining = set(diagnostics)
                restored = False
                # The restore must occur inside the same explicitly approved method.
                for restore_line in lines[index + 1:(body_end - 1 if body_end else index + 1)]:
                    restore = re.match(r"\s*#pragma\s+warning\s+restore\s*(.*)", restore_line)
                    if restore:
                        restored_ids = diagnostic_ids(restore.group(1).split("//")[0])
                        if not restored_ids:
                            remaining.clear()
                        else:
                            remaining.difference_update(restored_ids)
                        if not remaining:
                            restored = True
                            break
                if not symbol or not restored:
                    errors.append(f"{relative}:{index + 1}: warning suppression extends outside a recognized method or lacks matching restore within approved symbol")
            found.append({"file": relative, "line": index + 1, "diagnostics": diagnostics,
                          "rule": marker.group(1) if marker else None, "exception": marker.group(2) if marker else None,
                          "symbol": symbol})

        if path.suffix == ".cs":
            for index, line in enumerate(lines):
                pragma = re.search(r"^\s*#pragma\s+warning\s+disable\s*(.*)", line)
                if pragma:
                    add(index, diagnostic_ids(pragma.group(1).split("//")[0]) or ["ALL"], body=True)
            # Detection spans lines: assembly/global SuppressMessage attributes must not evade lint.
            for match in re.finditer(r"\[\s*(?:(?:assembly|module)\s*:\s*)?(?:[\w.]+\.)?(?:UnconditionalSuppressMessage|SuppressMessage)(?:Attribute)?\s*\(", text):
                add(text[:match.start()].count("\n"), ["ATTRIBUTE"])
        elif path.name in {".editorconfig", ".globalconfig"}:
            for index, line in enumerate(lines):
                if re.search(r"^\s*dotnet_(?:diagnostic|analyzer_diagnostic)\..*severity\s*=\s*(?:none|silent)", line):
                    add(index, ["EDITORCONFIG"])
        else:
            # MSBuild XML values may span lines or use imports/property expansions.
            tree = ET.fromstring(text)
            for element in tree.iter():
                if element.tag.rsplit("}", 1)[-1] not in {"NoWarn", "WarningsNotAsErrors"}:
                    continue
                value = "".join(element.itertext())
                value = re.sub(r"\$\((NoWarn|WarningsNotAsErrors)\)", "", value)
                if not value.strip(" \n\t\r;"):
                    continue
                tag = element.tag.rsplit("}", 1)[-1]
                index = next((index for index, line in enumerate(lines) if re.search(r"<" + tag + r"\b", line)), 0)
                add(index, diagnostic_ids(value) or ["CONFIGURATION"])
        for index, line in enumerate(lines):
            if "OVERRIDE(" in line and not MARKER.search(line):
                errors.append(f"{relative}:{index + 1}: malformed OVERRIDE reference")
    return found, errors


def lint(root: Path, now: datetime, evidence: dict | None = None) -> list[str]:
    catalog, errors = validate_catalog(root)
    registry_path = root / "docs/standards/overrides.yaml"
    registry = load_json(registry_path)
    digest = sha256(registry_path.read_bytes())
    if registry.get("schema_version") != 1 or registry.get("standard_version") != STANDARD_VERSION:
        errors.append("DOTNET-GOV-004: unsupported override schema/version")
    records = registry.get("overrides", [])
    ids = [record.get("id") for record in records]
    if len(set(ids)) != len(ids):
        errors.append("DOTNET-GOV-004: duplicate exception IDs")
    current, scan_errors = findings(root)
    errors.extend(scan_errors)
    for path in source_paths(root):
        for marker in MARKER.finditer(path.read_text(encoding="utf-8")):
            if marker.group(1) not in {item["id"] for item in catalog["rules"]} or marker.group(2) not in ids:
                errors.append(f"{path.relative_to(root)}: unknown OVERRIDE rule/exception {marker.group(0)}")
    used = set()
    for finding in current:
        record = next((item for item in records if item.get("id") == finding["exception"]), None)
        if record is None:
            errors.append(f"DOTNET-REPO-004 {finding['file']}:{finding['line']}: unregistered suppression {finding['diagnostics']} (independent approval required)")
            continue
        used.add(record["id"])
        errors.extend(validate_record(record, catalog, now, digest, evidence))
        if finding["rule"] != record.get("rule") or finding["file"] not in record.get("scope", {}).get("files", []) or not set(finding["diagnostics"]).issubset(record.get("diagnostics", [])) or finding["symbol"] not in record.get("scope", {}).get("symbols", []):
            errors.append(f"{record['id']} {record['rule']} {finding['file']}:{finding['line']}: marker/diagnostic/file/symbol scope mismatch")
    # Registry records must also be well formed when currently unused; active unused records are errors until retired.
    for record in records:
        if record.get("id") not in used:
            if record.get("status") == "active":
                errors.append(f"{record.get('id')} {record.get('rule')}: unused active record must be retired")
            elif record.get("status") in {"pending", "retired"}:
                # Validate structural/rule/scope/date fields even for unused proposals/history.
                structural = validate_record(record, catalog, now, digest, None)
                ignored = ("record cannot authorize", "missing independent approval", "missing HTTPS remediation", "missing protected external", "expired at")
                errors.extend(error for error in structural if not any(text in error for text in ignored))
            else:
                errors.append(f"{record.get('id')}: unknown exception status")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--release", action="store_true", help="release/trusted job; evidence cannot come from repository")
    parser.add_argument("--trusted-evidence", type=Path)
    parser.add_argument("--evidence-sha256", help="digest supplied by independently protected job, never PR input")
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        evidence = None
        if args.trusted_evidence:
            path = args.trusted_evidence.resolve()
            if not args.release or path.is_relative_to(root) or not args.evidence_sha256 or sha256(path.read_bytes()) != args.evidence_sha256:
                raise ValueError("external evidence requires release mode, outside-repository path and protected expected digest")
            evidence = load_json(path)
        errors = lint(root, datetime.now(timezone.utc), evidence)
    except (ValueError, TypeError, KeyError, AttributeError, OSError, ET.ParseError) as error:
        errors = [f"DOTNET-GOV-001/004: malformed/missing input: {error}"]
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print("Adopted 43-rule catalog and local override lint passed; external qualification/review are separate controls.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
