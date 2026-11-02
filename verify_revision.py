#!/usr/bin/env python3
"""Reproduce the additional exact checks for the revised manuscript.

First generate tree_optima.csv with optimize_trees.cpp. The independent
small-order optimizer below uses full eight-coordinate polynomial states;
the C++ program uses parity-compressed four-coordinate states.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction
import hashlib
import itertools as it
import json
import math
from pathlib import Path
import platform
import time
from verify import (P, exact_F, leaf_count, matching_number, reduced_count,
                    star_four, write_csv)


def parse_tree(word):
    adj, stack = [], []
    for ch in word:
        if ch == '(':
            v = len(adj)
            adj.append(set())
            if stack:
                adj[v].add(stack[-1])
                adj[stack[-1]].add(v)
            stack.append(v)
        elif ch == ')':
            assert stack
            stack.pop()
        else:
            raise ValueError('Invalid tree encoding')
    assert not stack and len(adj) > 1
    assert sum(map(len, adj)) == 2*(len(adj)-1)
    return adj


def python_optima(max_n=31):
    def mul(a, b):
        return tuple(sum(a[i]*b[k-i] for i in range(k+1)) for k in range(4))
    def combine(f, t):
        a, b, c, d = f[:4], f[4:], t[:4], t[4:]
        return mul(a, d) + tuple(x+y for x, y in zip(mul(b, d), mul(a, c)))
    def prune(values):
        kept = []
        for v in sorted(values, reverse=True):
            if not any(all(x >= y for x, y in zip(w, v)) for w in kept):
                kept.append(v)
        return kept
    forests = [[(1,0,0,0,0,0,0,0)]]
    trees = [[], [(1,0,0,0,0,1,0,0)]]
    maxima = {2: 0}
    for m in range(1, max_n-1):
        states = {combine(f,t) for j in range(1,m+1)
                  for f in forests[m-j] for t in trees[j]}
        forests.append(prune(states))
        trees.append(forests[-1])
        n = m+2
        maxima[n] = max(t[3]+t[6] for t in trees[n-1])
    return maxima


def spider_value(n):
    R = (n-5)//2
    if n % 2 == 0 or R < 4:
        return None
    q, r = divmod(R,4)
    return 4*q**3 + 3*r*q*q + r*(r-1)*q + r*(r-1)*(r-2)//6


def compact_example():
    # Two degree-three branch vertices. The three weight-one outer slots
    # are bare leaves; the middle slot has weight 2 and the fourth outer 3.
    adj = [set() for _ in range(7)]
    def edge(u,v):
        adj[u].add(v); adj[v].add(u)
    # 0,1 branches; 2 middle; 3,4,5 bare leaves; 6 outer centre.
    for u,v in [(0,2),(1,2),(0,3),(0,4),(1,5),(1,6)]: edge(u,v)
    for centre, arms in [(2,2),(6,3)]:
        for _ in range(arms):
            v=len(adj); adj.extend([set(),set()])
            edge(centre,v); edge(v,v+1)
    return adj


def run(optima, output):
    started=time.time()
    output.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader(optima.open()))
    values={int(r['n']):int(r['max_s3']) for r in rows}
    assert sorted(values)==list(range(2,max(values)+1))
    witness_checks=0
    comparison=[]
    for row in rows:
        n=int(row['n']); value=int(row['max_s3'])
        adj=parse_tree(row['witness_parentheses'])
        assert len(adj)==n and leaf_count(adj)==value
        deficit=n-2*matching_number(adj)
        if deficit==3:
            assert reduced_count(adj,{v for v,ns in enumerate(adj) if len(ns)==1})[0]==value
        witness_checks+=1
        if n % 2:
            S=(n-7)//2
            f, hs=exact_F(S) if S>=0 else (None,[])
            comparison.append(dict(n=n,M3=value,F=f,
                four_spider=spider_value(n),optimal_h=';'.join(map(str,hs)),
                witness_deficiency=deficit))
    reference=python_optima(min(31,max(values)))
    assert all(values[n]==v for n,v in reference.items())
    base=Path(__file__).parent/'results/tree_enumeration.csv'
    for row in csv.DictReader(base.open()):
        n=int(row['n'])
        if n in values: assert values[n]==int(row['max_s3'])
    # All new allocations with totals 41..60; the original program checks 0..40.
    allocation_checks=0
    for S in range(41,61):
        brute=0
        for cuts in it.combinations(range(S+4),4):
            a,b,h,c,d=(cuts[0],cuts[1]-cuts[0]-1,cuts[2]-cuts[1]-1,
                       cuts[3]-cuts[2]-1,S+3-cuts[3])
            val=P(a,b,h,c,d)
            assert 243*val<=16*S**3
            brute=max(brute,val)
            allocation_checks+=1
        assert brute==exact_F(S)[0]
    # Table 1: for each parity, both sides are polynomials of degree at most
    # 3 in each of r,s and at most 1 in h. This interpolation grid certifies
    # the identities exactly, beyond merely sampling the allowed region.
    parity_checks=0
    def g(h,x,z): return z*(x*x//4)+x*(z*z//4)+h*x*z
    for e,k in it.product(range(2),repeat=2):
        for r,s,h in it.product(range(4),range(4),range(2)):
            x,z=2*r+e,2*s+k
            table={(0,0):(r+s)*(r-s-1)+h*(2*(r-s)-1),
                   (1,1):(r-s)*(r+s+1)+h*(2*(r-s)-1),
                   (0,1):(r-s-1)*(r+s+1)+2*h*(r-s-1),
                   (1,0):(r-s)*(r+s)+2*h*(r-s)}
            assert g(h,x-1,z+1)-g(h,x,z)==table[e,k]
            parity_checks+=1
    poly=[-265599,24543,-12501,13]
    def numerator(s):return sum(c*s**i for i,c in enumerate(poly))
    for S in range(4):
        lhs=Fraction(16*(S-9)**3,243)-Fraction((S+5)**3,16)-Fraction((S+5)**2,2)
        assert lhs==Fraction(numerator(S),3888)
    # All coefficients after S=1000+u are positive: positivity for every u>=0.
    shifted=[sum(poly[i]*math.comb(i,j)*1000**(i-j) for i in range(j,4)) for j in range(4)]
    assert all(c>0 for c in shifted)
    assert numerator(959)<0<numerator(960)
    # Integer four-spider formula at every feasible odd order up to 2007.
    spider_checks=0
    for n in range(13,2008,2):
        q,r=divmod((n-5)//2,4)
        arms=[q+(i<r) for i in range(4)]
        direct=sum(math.prod(t) for t in it.combinations(arms,3))
        assert direct==spider_value(n)
        if n<=151:assert leaf_count(star_four(arms))==direct
        spider_checks+=1
    example=compact_example()
    assert len(example)==17 and leaf_count(example)==26
    assert len(example)-2*matching_number(example)==3
    # Counts, not the classification of every maximizer, are checked here.
    threshold_test=[n for n in values if n>=121 and n%2]
    if 119 in values:
        assert values[119]>exact_F(56)[0]
    assert all(values[n]==exact_F((n-7)//2)[0] for n in threshold_test)
    write_csv(output/'small_order_comparison.csv',comparison)
    summary=dict(python=platform.python_version(),max_n=max(values),
        witness_checks=witness_checks,python_reference_max_n=max(reference),
        python_reference_orders=len(reference),additional_allocation_checks=allocation_checks,
        total_allocations_through_60=math.comb(65,5),parity_identity_grid_checks=parity_checks,
        separation_polynomial_shift_1000=shifted,spider_integer_checks=spider_checks,
        equality_with_F_verified_odd_range=[121,max(threshold_test)] if threshold_test else None,
        classification_verified_by_optimizer=False,failures=0,
        elapsed_seconds=round(time.time()-started,3),
        source_sha256={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                       for name in ['verify.py','optimize_trees.cpp','verify_revision.py']})
    (output/'revision_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--optima',type=Path,default=Path('results/tree_optima.csv'))
    parser.add_argument('--output',type=Path,default=Path('run-results'))
    args=parser.parse_args()
    run(args.optima,args.output)
