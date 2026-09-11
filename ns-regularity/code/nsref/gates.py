"""Evaluate a run against the preregistered gates in prereg/gates.json."""
import json, os, math

_OPS = {"<=": lambda a, b: a <= b, ">=": lambda a, b: a >= b,
        "<": lambda a, b: a < b, ">": lambda a, b: a > b,
        "==": lambda a, b: a == b}

DEFAULT_GATES = os.path.join(os.path.dirname(__file__), "..", "..", "prereg", "gates.json")


def load_gates(path=None):
    with open(os.path.abspath(path or DEFAULT_GATES)) as f:
        return json.load(f)


def evaluate(measured, gates_doc=None, applicable=None):
    """measured: {quantity_name: value}. Returns per-gate PASS/FAIL/NOT_MEASURED.

    A quantity that was not measured is NOT a pass. Gates fail closed.
    """
    doc = gates_doc or load_gates()
    results = {}
    for name, spec in doc["gates"].items():
        if applicable is not None and name not in applicable:
            results[name] = {"status": "NOT_APPLICABLE", "node": spec["node"], "criteria": []}
            continue
        crits, status = [], "PASS"
        for c in spec["criteria"]:
            q, op, want = c["quantity"], c["op"], c["value"]
            if q not in measured:
                crits.append({**c, "measured": None, "result": "NOT_MEASURED"})
                status = "FAIL" if status != "FAIL" else status
                continue
            got = measured[q]
            if isinstance(got, float) and math.isnan(got):
                ok = False
            else:
                try:
                    ok = _OPS[op](got, want)
                except TypeError:
                    ok = False
            crits.append({**c, "measured": got, "result": "PASS" if ok else "FAIL"})
            if not ok:
                status = "FAIL"
        results[name] = {"status": status, "node": spec["node"], "criteria": crits}
    n_pass = sum(1 for r in results.values() if r["status"] == "PASS")
    n_appl = sum(1 for r in results.values() if r["status"] != "NOT_APPLICABLE")
    return {"gates": results, "n_pass": n_pass, "n_applicable": n_appl,
            "all_pass": n_pass == n_appl and n_appl > 0}


def format_table(res):
    lines = [f"{'GATE':<34} {'NODE':<8} STATUS", "-" * 56]
    for name, r in res["gates"].items():
        lines.append(f"{name:<34} {r['node']:<8} {r['status']}")
        for c in r["criteria"]:
            if c["result"] != "PASS":
                m = c["measured"]
                ms = "not measured" if m is None else (f"{m:.4g}" if isinstance(m, float) else str(m))
                lines.append(f"    - {c['quantity']} {c['op']} {c['value']}  (got {ms})")
    lines.append("-" * 56)
    lines.append(f"{res['n_pass']}/{res['n_applicable']} applicable gates passed"
                 f"   ALL_PASS={res['all_pass']}")
    return "\n".join(lines)
