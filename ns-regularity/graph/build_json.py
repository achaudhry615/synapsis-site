"""Regenerate proof_obligations.json from nodes.csv + edges.csv (single source of truth)."""
import csv, json, collections, hashlib, pathlib

here = pathlib.Path(__file__).parent
nodes = list(csv.DictReader(open(here / "nodes.csv")))
edges = list(csv.DictReader(open(here / "edges.csv")))

ids = {n["id"] for n in nodes}
dangling = [(e["src"], e["dst"]) for e in edges if e["src"] not in ids or e["dst"] not in ids]
if dangling:
    raise SystemExit(f"dangling edges: {dangling}")

out_edges = collections.defaultdict(list)
in_edges = collections.defaultdict(list)
for e in edges:
    out_edges[e["src"]].append({"to": e["dst"], "relation": e["relation"], "note": e["note"]})
    in_edges[e["dst"]].append({"from": e["src"], "relation": e["relation"], "note": e["note"]})

doc = {
    "schema": "ns-regularity/proof-obligation-graph/v1",
    "disclaimer": (
        "3D Navier-Stokes global regularity is OPEN. No node in this graph is a proof of it. "
        "Status values are descriptive bookkeeping, not claims of resolution."
    ),
    "status_vocabulary": {
        "CERTIFIED_STANDARD": "published mathematics, cited, not our contribution",
        "CERTIFIED_STANDARD_VERIFY_CITATION": "as above but the exact hypotheses/constants were not verifiable from this runtime",
        "PROVED_ELEMENTARY": "proved here; elementary",
        "NO_GO_PROVEN": "proved here that this route cannot close",
        "THEOREM_CONDITIONAL": "proved here by reduction; hypothesis unproved",
        "CONJECTURE": "stated so it can be false; not assumed anywhere",
        "OPEN": "no path proposed",
        "DIAGNOSTIC_ONLY": "measurable, but no theorem connects it to regularity",
        "NUMERICAL_GATE": "preregistered pass/fail criterion",
        "PREREGISTERED": "falsifiable hypothesis fixed before the campaign",
        "BLOCKED": "cannot proceed until a prerequisite changes status",
        "SMOKE_TEST_ONLY": "pipeline validation; carries no evidential weight",
    },
    "nodes": [dict(n, out_edges=out_edges[n["id"]], in_edges=in_edges[n["id"]]) for n in nodes],
    "summary": dict(collections.Counter(n["status"] for n in nodes)),
    "open_problems": [n["id"] for n in nodes if n["status"] == "OPEN"],
    "conditional_theorems": [n["id"] for n in nodes if n["status"] == "THEOREM_CONDITIONAL"],
}
text = json.dumps(doc, indent=2)
(here / "proof_obligations.json").write_text(text + "\n")
print("nodes", len(nodes), "edges", len(edges))
print("sha256", hashlib.sha256(text.encode()).hexdigest()[:16])
print("open:", doc["open_problems"])
print("conditional:", doc["conditional_theorems"])
