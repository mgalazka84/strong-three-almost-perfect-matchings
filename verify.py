#!/usr/bin/env python3
"""Exact, reproducible checks for strong three-almost-perfect matchings.

No floating point arithmetic is used in the verification. NetworkX is used
only to generate nonisomorphic trees and as an independent matching solver.
The counting and Gallai--Edmonds reduction routines are implemented here.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools as it
import json
import math
from pathlib import Path
import random
import sys
import time
import networkx as nx


def adjacency(g):
    return [set(g.neighbors(v)) for v in range(len(g))]


def matching_number(adj, removed=frozenset()):
    """Linear-time maximum matching size in a forest, by rooted DP."""
    seen = set(removed)
    answer = 0
    for root in range(len(adj)):
        if root in seen:
            continue
        parent = {root: -1}
        order = [root]
        seen.add(root)
        for v in order:
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    parent[w] = v
                    order.append(w)
        free, closed = {}, {}
        for v in reversed(order):
            children = [w for w in adj[v] if parent.get(w) == v]
            free[v] = sum(closed[w] for w in children)
            closed[v] = free[v] + max(
                [0] + [1 + free[w] - closed[w] for w in children]
            )
        answer += closed[root]
    return answer


def coefficients(adj, eligible, degree=3):
    """Coefficient j counts matchings with j monomers, all in eligible."""
    def mul(a, b):
        c = [0] * (degree + 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b[:degree + 1 - i]):
                    c[i+j] += x*y
        return c
    one = [1] + [0]*degree
    zero = [0]*(degree+1)
    seen = set()
    answer = one
    for root in range(len(adj)):
        if root in seen:
            continue
        parent = {root: -1}
        order = [root]
        seen.add(root)
        for v in order:
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    parent[w] = v
                    order.append(w)
        reserved, closed = {}, {}
        for v in reversed(order):
            free, paired = one[:], zero[:]
            for w in adj[v]:
                if parent.get(w) != v:
                    continue
                paired = [a+b for a, b in zip(
                    mul(paired, closed[w]), mul(free, reserved[w]))]
                free = mul(free, closed[w])
            reserved[v] = free
            closed[v] = paired
            if v in eligible:
                for j in range(degree):
                    closed[v][j+1] += free[j]
        answer = mul(answer, closed[root])
    return answer


def leaf_count(adj):
    return coefficients(adj, {v for v, ns in enumerate(adj) if len(ns) == 1})[3]


def direct_leaf_count(g):
    """Independent: delete each leaf triple and use general matching."""
    leaves = [v for v in g if g.degree(v) == 1]
    result = 0
    for removed in it.combinations(leaves, 3):
        h = g.copy()
        h.remove_nodes_from(removed)
        m = nx.max_weight_matching(h, maxcardinality=True)
        result += 2*len(m) == len(h)
    return result


def reduced_count(adj, eligible):
    """Count eligible maximum matchings through the weighted reduced forest.

    D is obtained from independently computed vertex-deletion matching sizes.
    The polynomial is evaluated by enumerating branch-to-slot matchings,
    rather than by using any specialized deficiency-three formula.
    """
    n = len(adj)
    nu = matching_number(adj)
    D = {v for v in range(n) if matching_number(adj, {v}) == nu}
    A = set().union(*(adj[v] for v in D)) - D if D else set()
    C = set(range(n)) - D - A
    assert all(not (adj[v] & D) for v in D)
    assert all(len(adj[v] & D) >= 2 for v in A)
    assert n - 2*nu == len(D)-len(A)
    assert matching_number(adj, D | A)*2 == len(C)
    B = {v for v in A if len(adj[v] & D) >= 3}
    keep = D | (A-B)
    seen = set()
    slots, vertex_slot = [], {}
    for root in sorted(keep):
        if root in seen:
            continue
        comp = {root}
        todo = [root]
        seen.add(root)
        for v in todo:
            for w in adj[v] & keep:
                if (v in D) == (w in D):
                    continue
                if w not in seen:
                    seen.add(w)
                    comp.add(w)
                    todo.append(w)
        idx = len(slots)
        ds = comp & D
        assert len(ds) == len(comp - D) + 1
        slots.append(len(ds & eligible))
        for v in ds:
            vertex_slot[v] = idx
    branch_neighbors = []
    for v in sorted(B):
        nei = [vertex_slot[w] for w in adj[v] & D]
        assert len(nei) == len(set(nei))
        branch_neighbors.append(nei)
    total = 0
    def rec(i, used):
        nonlocal total
        if i == len(branch_neighbors):
            total += math.prod(x for j, x in enumerate(slots) if j not in used)
            return
        for j in branch_neighbors[i]:
            if j not in used:
                rec(i+1, used | {j})
    rec(0, set())
    d = n-2*nu
    if d:
        assert len(B)+len(slots) <= 3*d-2
    return total, d, len(eligible & D), len(B), len(slots)


def P(a, b, y, c, d):
    return a*b*(c+d) + c*d*(a+b) + y*(a+b)*(c+d)


def H(S, y):
    x = (S-y)//2
    z = S-y-x
    return z*(x*x//4) + x*(z*z//4) + y*x*z


def exact_F(S):
    values = [H(S, y) for y in range(S+1)]
    best = max(values)
    return best, [y for y, v in enumerate(values) if v == best]


def decorated(arms, edges):
    adj = [set() for _ in arms]
    def edge(u, v):
        adj[u].add(v)
        adj[v].add(u)
    for u, v in edges:
        edge(u, v)
    for v, a in enumerate(arms):
        for _ in range(a):
            w = len(adj)
            adj.extend([set(), set()])
            edge(v, w)
            edge(w, w+1)
    return adj


def double_fork(a, b, y, c, d):
    return decorated([0, 0, a, b, y, c, d],
                     [(0,2), (0,3), (0,4), (1,4), (1,5), (1,6)])


def star_four(arms):
    return decorated([0]+list(arms), [(0,i) for i in range(1,5)])


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def run(output, max_n=17, kernel_n=12, direct_n=11):
    output.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261030)
    started = time.time()
    rows = []
    kernel_checks = direct_checks = trees_checked = 0
    for n in range(2, max_n+1):
        best = -1
        best_count = count = 0
        by_def = {}
        for g in nx.nonisomorphic_trees(n):
            count += 1
            trees_checked += 1
            adj = adjacency(g)
            leaves = {v for v in g if g.degree(v) == 1}
            s3 = leaf_count(adj)
            d = n - 2*matching_number(adj)
            by_def[d] = max(by_def.get(d,0), s3)
            if n % 2:
                assert 243*s3 <= 2*(n+3)**3
            else:
                assert s3 == 0
            if d == 1:
                assert 162*s3 <= (n+1)**3
            if s3 > best:
                best, best_count = s3, 1
            elif s3 == best:
                best_count += 1
            if n <= direct_n:
                assert direct_leaf_count(g) == s3
                direct_checks += 1
            if n <= kernel_n:
                subsets = [set(g), leaves, set()]
                if n <= 7:
                    subsets = [{v for v in g if mask>>v & 1}
                               for mask in range(1<<n)]
                else:
                    subsets += [{v for v in g if rng.randrange(2)} for _ in range(5)]
                for allowed in subsets:
                    value, deficit, L, nb, ns = reduced_count(adj, allowed)
                    dp = coefficients(adj, allowed, deficit)[deficit]
                    assert dp == value
                    if deficit == 1:
                        assert value <= L
                    if deficit == 2:
                        assert 3*value <= L*L
                    if deficit == 3:
                        assert 243*value <= 16*L**3
                    kernel_checks += 1
        rows.append(dict(n=n, unlabelled_trees=count, max_s3=best,
                         maximizers=best_count,
                         max_deficiency_1=by_def.get(1,0),
                         max_deficiency_3=by_def.get(3,0)))
        print(json.dumps(rows[-1]), flush=True)
    write_csv(output/'tree_enumeration.csv', rows)
    # All five-variable integer allocations, compared with the one-dimensional rule.
    allocation_checks = 0
    for S in range(41):
        brute = 0
        for cuts in it.combinations(range(S+4), 4):
            a = cuts[0]
            b = cuts[1]-cuts[0]-1
            y = cuts[2]-cuts[1]-1
            c = cuts[3]-cuts[2]-1
            d = S+3-cuts[3]
            val = P(a,b,y,c,d)
            assert 243*val <= 16*S**3
            if min(a,b,y,c,d) <= 1:
                assert 16*val <= S**3+8*S*S
            brute = max(brute, val)
            allocation_checks += 1
        assert brute == exact_F(S)[0]
    # Independently evaluate the proposed family with the polynomial DP.
    family_rows = []
    for t in [1,2,3,4,5,6,7,8,10,20,50,112,200]:
        adj = double_fork(2*t,2*t,t,2*t,2*t)
        value = leaf_count(adj)
        assert value == 48*t**3
        n = len(adj)
        S_old = (n-5)//2
        q,r = divmod(S_old,4)
        arms = [q+(i<r) for i in range(4)]
        old = leaf_count(star_four(arms))
        assert old == sum(math.prod(c) for c in it.combinations(arms,3))
        assert exact_F(9*t)[0] == value
        if t >= 7:
            assert value > old
        family_rows.append(dict(t=t,n=n,double_fork=value,
                                balanced_four_spider=old,difference=value-old))
    write_csv(output/'explicit_family.csv', family_rows)
    # Finite formula values; not an assertion about all small trees.
    formula_rows=[]
    rounded_formula_instances = 0
    for S in list(range(4,101))+[1000,1001,1002,1008,1500,2000]:
        value, ys = exact_F(S)
        formula_rows.append(dict(n=2*S+7,S=S,F=value,
                                 optimal_y=';'.join(map(str,ys)),
                                 global_theorem_applies=(S>=1000)))
        if S >= 1000:
            for y in ys:
                x = (S-y)//2
                z = S-y-x
                a, b = x//2, x-x//2
                c, d = z//2, z-z//2
                assert min(a,b,y,c,d) >= 2
                adj = double_fork(a,b,y,c,d)
                assert len(adj) == 2*S+7
                assert len(adj)-2*matching_number(adj) == 3
                assert leaf_count(adj) == value
                rounded_formula_instances += 1
    write_csv(output/'integer_formula.csv', formula_rows)
    # Certificate for the numerical separation in the eventual theorem.
    for S in range(1000,10001):
        assert 13*S**3-12501*S**2+24543*S-265599 > 0
    summary=dict(python=sys.version.split()[0], networkx=nx.__version__,
        max_n=max_n,kernel_n=kernel_n,direct_n=direct_n,
        unlabelled_trees=trees_checked,independent_triple_checks=direct_checks,
        weighted_reduction_checks=kernel_checks,
        integer_allocations_checked=allocation_checks,
        family_instances=len(family_rows),failures=0,
        rounded_formula_instances=rounded_formula_instances,
        elapsed_seconds=round(time.time()-started,3),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('results'))
    parser.add_argument('--max-n',type=int,default=17)
    parser.add_argument('--kernel-n',type=int,default=12)
    parser.add_argument('--direct-n',type=int,default=11)
    args=parser.parse_args()
    run(args.output,args.max_n,args.kernel_n,args.direct_n)
