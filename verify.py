#!/usr/bin/env python3
"""Independent checker: recompute Pg and |CC| from the JSON the backend returns (pure python, same formulas as the paper / frontend)."""
import json, math, sys, collections
d = json.load(sys.stdin); p = d["params"]; R, C = d["rows"], d["cols"]
S = [tuple(x) for x in d["sensors"]]; sink = tuple(d["sink"])
assert len(set(S)) == len(S) and sink not in S, "duplicate / sink overlap"
def psg(x):
    if x <= p["Rcon"]: return p["alpha"]
    if x <= p["Rmax"]: return 1 / x ** p["beta"]
    return 0
bad = 0; mn = 1
for r in range(1, R + 1):
    for c in range(1, C + 1):
        miss = 1
        for (a, b) in S:
            miss *= 1 - psg(math.hypot(r - a, c - b))
        mn = min(mn, 1 - miss)
        if 1 - miss < p["Tg"] - 1e-9: bad += 1
nodes = S + [sink]; seen = {sink}; q = collections.deque([sink])
while q:
    u = q.popleft()
    for v in nodes:
        if v not in seen and math.hypot(u[0]-v[0], u[1]-v[1]) <= p["Rt"]: seen.add(v); q.append(v)
conn = len(seen) - 1
print(f"sensors={len(S)} cells_below_Tg={bad} Pg_min={mn:.4f} connected={conn}/{len(S)} (claimed feasible={d['feasible']})")
ok = bad == 0 and (conn == len(S) or not p["conn"])
print("VERIFIED" if ok == d["feasible"] else "MISMATCH", "-", "constraints satisfied" if ok else "constraints violated")
