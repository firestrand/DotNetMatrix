"""Tests for maintainer approval, exact exception scope and standards lint."""
from __future__ import annotations

import copy
import importlib.util
import json
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("standards_gate", Path(__file__).resolve().parents[1] / "check-standards.py")
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)
ROOT = Path(__file__).resolve().parents[2]


class StandardsGateTests(unittest.TestCase):
    def setUp(self):
        self.catalog = GATE.load_json(ROOT / "docs/standards/dotnet-rules.json")
        self.record = copy.deepcopy(GATE.load_json(ROOT / "docs/standards/overrides.yaml")["overrides"][0])
        self.record.update(status="active", approved_by=self.record["owner"])
        self.now = datetime(2026, 10, 3, tzinfo=timezone.utc)

    def validate(self):
        return GATE.validate_record(self.record, self.catalog, self.now)

    def test_exact_source_and_43_catalog_ids(self):
        catalog, errors = GATE.validate_catalog(ROOT)
        self.assertEqual([], errors)
        self.assertEqual(43, len(catalog["rules"]))
        group = next(rule for rule in catalog["rules"] if rule["id"] == "DOTNET-BASE-004")
        self.assertGreater(group["source"]["end_line"], 122)

    def test_maintainer_can_approve_without_external_evidence(self):
        self.assertEqual(self.record["owner"], self.record["approved_by"])
        self.assertNotIn("remediation_issue", self.record)
        self.assertEqual([], self.validate())

    def test_utc_expiration_is_exclusive_at_midnight_when_specified(self):
        self.record["expires_on"] = "2026-11-01"
        self.now = datetime(2026, 11, 1, tzinfo=timezone.utc)
        self.assertTrue(any("expired" in error for error in self.validate()))
        self.now = datetime(2026, 10, 31, 23, 59, 59, tzinfo=timezone.utc)
        self.assertEqual([], self.validate())

    def test_leap_date_and_malformed_date(self):
        self.record.update(created_on="2028-02-28", expires_on="2028-02-29")
        self.now = datetime(2028, 2, 28, 23, tzinfo=timezone.utc)
        self.assertEqual([], self.validate())
        self.record["expires_on"] = "2027-02-29"
        self.assertTrue(any("calendar" in error for error in self.validate()))

    def test_unknown_rule_clause_id_and_required_fields(self):
        for field, value, expected in [("rule", "DOTNET-FAKE-999", "unknown rule"),
                                       ("clause", "entire-rule-group", "clause selector"),
                                       ("id", "bad-id", "malformed exception")]:
            with self.subTest(field=field):
                saved = self.record[field]
                self.record[field] = value
                self.assertTrue(any(expected in error for error in self.validate()))
                self.record[field] = saved
        del self.record["security_impact"]
        self.assertTrue(any("missing registry fields" in error for error in self.validate()))

    def test_scope_cannot_be_wildcard_absolute_or_traversal(self):
        for path in ["**/*.cs", "/source/file.cs", "../file.cs"]:
            self.record["scope"]["files"] = [path]
            self.assertTrue(any("exact repository-relative" in error for error in self.validate()))

    def test_missing_approval_and_inactive_records_are_rejected(self):
        for approver in [None, "", "  "]:
            with self.subTest(approver=approver):
                self.record["approved_by"] = approver
                self.assertTrue(any("maintainer approval" in error for error in self.validate()))
        self.record["approved_by"] = self.record["owner"]
        for status in ["retired", "pending"]:
            self.record["status"] = status
            self.assertTrue(any(f"{status} record" in error for error in self.validate()))

    def test_permanent_approval_has_no_artificial_expiration(self):
        self.record["expires_on"] = None
        self.now = datetime(2040, 1, 1, tzinfo=timezone.utc)
        self.assertEqual([], self.validate())

    def test_maintainer_approval_does_not_authorize_expanded_diagnostics(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "docs/standards"
            target.mkdir(parents=True)
            for name in ["engineering-standard.md", "dotnet-rules.json"]:
                (target / name).write_bytes((ROOT / "docs/standards" / name).read_bytes())
            record = copy.deepcopy(self.record)
            record["scope"] = {"files": ["Test.cs"], "symbols": ["Example.Tests.Method"]}
            (target / "overrides.yaml").write_text(json.dumps({"schema_version": 1, "standard_version": GATE.STANDARD_VERSION, "overrides": [record]}))
            (root / "Test.cs").write_text('namespace Example;\npublic class Tests {\n public void Method() {\n// OVERRIDE(DOTNET-REPO-004, EXC-2026-001)\n#pragma warning disable SYSLIB0050, CA1416\n#pragma warning restore SYSLIB0050, CA1416\n }\n}\n')
            self.assertTrue(any("scope mismatch" in error for error in GATE.lint(root, self.now)))

    def test_catalog_duplicate_unknown_and_source_drift_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "docs/standards"
            target.mkdir(parents=True)
            raw = (ROOT / "docs/standards/engineering-standard.md").read_bytes()
            (target / "engineering-standard.md").write_bytes(raw)
            for mutation in ["duplicate", "unknown", "range"]:
                catalog = copy.deepcopy(self.catalog)
                if mutation == "duplicate":
                    catalog["rules"][1]["id"] = catalog["rules"][0]["id"]
                elif mutation == "unknown":
                    catalog["rules"][1]["id"] = "DOTNET-UNKNOWN-999"
                else:
                    catalog["rules"][0]["source"]["end_line"] += 1
                (target / "dotnet-rules.json").write_text(json.dumps(catalog))
                self.assertTrue(GATE.validate_catalog(root)[1])
            catalog = copy.deepcopy(self.catalog)
            modified = raw + b"\nUnapproved source revision\n"
            catalog["source_sha256"] = GATE.sha256(modified)
            (target / "engineering-standard.md").write_bytes(modified)
            (target / "dotnet-rules.json").write_text(json.dumps(catalog))
            self.assertTrue(any("SHA-256" in error for error in GATE.validate_catalog(root)[1]))

    def test_pending_proposal_never_approves_live_suppression(self):
        proposal = copy.deepcopy(self.record)
        proposal.update(status="pending", approved_by=None)
        errors = GATE.validate_record(proposal, self.catalog, self.now)
        self.assertTrue(any("pending record" in error for error in errors))
        self.assertTrue(any("missing maintainer approval" in error for error in errors))

    def test_multiline_xml_attribute_and_globalconfig_suppressions_are_detected(self):
        fixtures = {
            "Project.csproj": "<Project><PropertyGroup>\n<NoWarn>\nCA1416\n</NoWarn>\n</PropertyGroup></Project>",
            "Assembly.cs": '[assembly:\nSystem.Diagnostics.CodeAnalysis.SuppressMessage("Security", "CA5350")]\n',
            ".globalconfig": "is_global = true\ndotnet_diagnostic.CA1416.severity = none\n",
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, content in fixtures.items():
                (root / name).write_text(content)
            findings, errors = GATE.findings(root)
            self.assertEqual([], errors)
            self.assertEqual(set(fixtures), {finding["file"] for finding in findings})
            self.assertTrue(all(finding["exception"] is None for finding in findings))
            self.assertEqual(["CA1416"], next(item for item in findings if item["file"] == "Project.csproj")["diagnostics"])

    def test_pragma_extent_cannot_cross_method_boundary_or_be_outside_method(self):
        fixtures = [
            'namespace Example;\npublic class Tests {\n public void First() {\n#pragma warning disable SYSLIB0050\n }\n public void Second() {\n#pragma warning restore SYSLIB0050\n }\n}\n',
            'namespace Example;\npublic class Tests {\n public void First() { }\n#pragma warning disable SYSLIB0050\n#pragma warning restore SYSLIB0050\n}\n',
            'namespace Example;\npublic class Tests {\n public void First() {\n#pragma warning disable SYSLIB0050, CA1416\n#pragma warning restore SYSLIB0050\n }\n}\n',
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "Tests.cs"
            for content in fixtures:
                with self.subTest(content=content):
                    path.write_text(content)
                    self.assertTrue(any("outside a recognized method" in error for error in GATE.findings(Path(temporary))[1]))
            path.write_text('namespace Example;\npublic class Tests {\n public void First() {\n#pragma warning disable SYSLIB0050, CA1416\n string literal = "}"; /* } */\n#pragma warning restore SYSLIB0050\n#pragma warning restore CA1416\n }\n}\n')
            findings, errors = GATE.findings(Path(temporary))
            self.assertEqual([], errors)
            self.assertEqual("Example.Tests.First", findings[0]["symbol"])

    def test_clause_cannot_select_only_a_fragment_of_a_normative_sentence(self):
        self.assertTrue(GATE.exact_normative_sentence(GATE.WARNING_CLAUSE, "**Warning policy.** " + GATE.WARNING_CLAUSE + " Next sentence."))
        self.assertFalse(GATE.exact_normative_sentence("publication **MUST** fail on the configured compiler/analyzer warnings.", GATE.WARNING_CLAUSE))

    def test_catalog_cannot_invent_a_clause_or_borrow_another_rule_clause(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "docs/standards"
            target.mkdir(parents=True)
            (target / "engineering-standard.md").write_bytes((ROOT / "docs/standards/engineering-standard.md").read_bytes())
            for clause in ["invented-waiver", GATE.WARNING_CLAUSE]:
                catalog = copy.deepcopy(self.catalog)
                catalog["rules"][0]["clauses"] = [clause]
                (target / "dotnet-rules.json").write_text(json.dumps(catalog))
                self.assertTrue(any("exact normative source sentence" in error for error in GATE.validate_catalog(root)[1]))

    def test_duplicate_json_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "registry.yaml"
            path.write_text('{"overrides":[], "overrides":[]}')
            with self.assertRaisesRegex(ValueError, "duplicate"):
                GATE.load_json(path)

    def test_unregistered_and_unknown_marker_scope_findings(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "Test.cs"
            path.write_text('namespace Example;\npublic class Tests {\n public void Method() {\n#pragma warning disable SYSLIB0050\n#pragma warning restore SYSLIB0050\n }\n}\n')
            found, errors = GATE.findings(root)
            self.assertEqual([], errors)
            self.assertEqual("Example.Tests.Method", found[0]["symbol"])
            self.assertIsNone(found[0]["exception"])
            path.write_text('// OVERRIDE(DOTNET-UNKNOWN-999, malformed)\n#pragma warning disable\n')
            found, errors = GATE.findings(root)
            self.assertEqual(["ALL"], found[0]["diagnostics"])
            self.assertTrue(any("malformed OVERRIDE" in error for error in errors))

    def test_duplicate_registry_ids_and_symbol_scope_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "docs/standards"
            target.mkdir(parents=True)
            for name in ["engineering-standard.md", "dotnet-rules.json"]:
                (target / name).write_bytes((ROOT / "docs/standards" / name).read_bytes())
            record = copy.deepcopy(self.record)
            record["scope"] = {"files": ["Test.cs"], "symbols": ["Example.Tests.WrongMethod"]}
            (target / "overrides.yaml").write_text(json.dumps({"schema_version": 1, "standard_version": GATE.STANDARD_VERSION, "overrides": [record, record]}))
            (root / "Test.cs").write_text('namespace Example;\npublic class Tests {\n public void Method() {\n// OVERRIDE(DOTNET-REPO-004, EXC-2026-001)\n#pragma warning disable SYSLIB0050\n#pragma warning restore SYSLIB0050\n }\n}\n')
            errors = GATE.lint(root, self.now)
            self.assertTrue(any("duplicate exception IDs" in error for error in errors))
            self.assertTrue(any("symbol scope mismatch" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
