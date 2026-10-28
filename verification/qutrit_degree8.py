"""Exact Drensky--Kasparian degree-eight central polynomial for 3x3 matrices.

Source: V. Drensky and A. Kasparian, Communications in Algebra 13(3),
745--752 (1985), DOI 10.1080/00927878508823188, theorem pp.746--747.

Words are mathematical matrix-product order, LEFT TO RIGHT. Python variable
labels are 0,...,7; exported labels are 1,...,8. A physical temporal order is
the reverse of the matrix-product word. No third-party dependencies are used.
"""
from collections import Counter
from itertools import permutations
from pathlib import Path
import argparse
import json
import random


def sign(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p))
                    for j in range(i+1, len(p))) % 2 else 1


def mul(a, b):
    out = Counter()
    for wa, ca in a.items():
        for wb, cb in b.items():
            out[wa+wb] += ca*cb
    return Counter({w:c for w,c in out.items() if c})


def add(*terms):
    out = Counter()
    for scale, poly in terms:
        for w,c in poly.items():
            out[w] += scale*c
    return Counter({w:c for w,c in out.items() if c})


def power(poly, n):
    result = Counter({():1})
    for _ in range(n):
        result = mul(result,poly)
    return result


def generate_coefficients_capelli():
    """The compact 20+5 placement formula on printed page747."""
    identity = Counter({():1})
    A = Counter({(4,5):1, (5,4):-1})
    B = Counter({(6,7):1, (7,6):-1})
    AB_BA = add((1,mul(A,B)), (1,mul(B,A)))
    out = Counter()

    def add_capelli(slots, scale):
        for p in permutations(range(4)):
            term = slots[0]
            for k in range(4):
                term = mul(mul(term, Counter({(p[k],):1})), slots[k+1])
            for w,c in term.items():
                out[w] += scale*sign(p)*c

    for i in range(5):
        for j in range(5):
            if i == j:
                continue
            slots = [identity]*5
            slots[i], slots[j] = A, B
            m = -1 if (i,j) in ((1,3),(3,1)) else 1
            add_capelli(slots,m)
    for k in range(5):
        slots = [identity]*5
        slots[k] = AB_BA
        n = -1 if k in (1,3) else 1
        add_capelli(slots,n)
    return {w:c for w,c in out.items() if c}


def generate_coefficients_exterior():
    """Independent expansion of the original Grassmann-algebra definition.

    F=U^4-XU^3X+U^2V+UVU+VU^2, U=X^2+Y^2+Z^2,
    V=X(Y^2+Z^2)X. The top exterior coefficient assigns each of
    x1,...,x4 to the four X positions and each of x5,x6 / x7,x8 to
    the two Y / Z positions; the Grassmann sign is the word's parity.
    """
    X = Counter({(0,):1})
    U = Counter({(0,0):1,(1,1):1,(2,2):1})
    V = Counter({(0,1,1,0):1,(0,2,2,0):1})
    F = add((1,power(U,4)),(-1,mul(mul(X,power(U,3)),X)),
            (1,mul(power(U,2),V)),(1,mul(mul(U,V),U)),
            (1,mul(V,power(U,2))))
    out = Counter()
    for pattern,c in F.items():
        if tuple(pattern.count(i) for i in range(3)) != (4,2,2):
            continue
        for px in permutations(range(4)):
            for py in permutations((4,5)):
                for pz in permutations((6,7)):
                    groups = [iter(px),iter(py),iter(pz)]
                    word = tuple(next(groups[g]) for g in pattern)
                    out[word] += c*sign(word)
    return {w:c for w,c in out.items() if c}


def mm(a,b):
    d=len(a)
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(d))
                       for j in range(d)) for i in range(d))


def madd(a,b):
    return tuple(tuple(x+y for x,y in zip(ra,rb)) for ra,rb in zip(a,b))


def mscale(c,a):
    return tuple(tuple(c*x for x in row) for row in a)


def ident(d):
    return tuple(tuple(int(i==j) for j in range(d)) for i in range(d))


def comm(a,b):
    return madd(mm(a,b),mscale(-1,mm(b,a)))


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def evaluate(coefficients, matrices):
    d=len(matrices[0])
    result=mscale(0,ident(d))
    for word,c in coefficients.items():
        value=ident(d)
        for label in word:
            value=mm(value,matrices[label])
        result=madd(result,mscale(c,value))
    return result


def scalar_trace_formula(matrices):
    """Scalar identity from printed page752; applies in dimension3."""
    A=comm(matrices[4],matrices[5])
    B=comm(matrices[6],matrices[7])
    return 2*sum(sign(p)*trace(mm(mm(matrices[p[0]],matrices[p[1]]),A))
                 *trace(mm(mm(matrices[p[2]],matrices[p[3]]),B))
                 for p in permutations(range(4)))


def unit(i,j,d=3):
    return tuple(tuple(int((a,b)==(i,j)) for b in range(d)) for a in range(d))


def source_example():
    """One-based source tuple: E11,E12,E23,E31,E12,E23,E32,E21."""
    return [unit(i,j) for i,j in ((0,0),(0,1),(1,2),(2,0),
                                 (0,1),(1,2),(2,1),(1,0))]


def verify(random_samples=12):
    coefficients=generate_coefficients_capelli()
    assert coefficients==generate_coefficients_exterior()
    assert len(coefficients)==2880
    assert Counter(coefficients.values())==Counter({-1:1440,1:1440})
    assert sum(abs(c) for c in coefficients.values())==2880
    assert all(sorted(w)==list(range(8)) for w in coefficients)
    assert evaluate(coefficients,source_example())==mscale(-2,ident(3))
    assert scalar_trace_formula(source_example())==-2
    rng=random.Random(20261028)
    scalars=[]
    for _ in range(random_samples):
        matrices=[tuple(tuple(rng.randrange(-2,3) for j in range(3))
                        for i in range(3)) for a in range(8)]
        scalar=scalar_trace_formula(matrices)
        assert evaluate(coefficients,matrices)==mscale(scalar,ident(3))
        scalars.append(scalar)
    # Lower-dimensional vanishing follows by embedding M2 as a corner of M3.
    matrices2=[tuple(tuple(rng.randrange(-2,3) for j in range(2))
                     for i in range(2)) for a in range(8)]
    assert evaluate(coefficients,matrices2)==mscale(0,ident(2))
    return dict(support_size=len(coefficients),l1_norm=2880,
                coefficient_histogram={str(k):v for k,v in
                                       sorted(Counter(coefficients.values()).items())},
                independent_generators_agree=True,
                source_matrix_unit_example_scalar=-2,
                exact_random_3x3_checks=random_samples,
                exact_random_scalars=scalars,
                one_exact_2x2_zero_check=True)


def export(path,coefficients):
    data=dict(name="Drensky--Kasparian qutrit central polynomial",
              source_doi="10.1080/00927878508823188",
              source_pages="746--747, 751--752",
              word_order="Mathematical matrix product, left to right; temporal order is reversed",
              labels="One-based x_1,...,x_8",
              degree=8,dimension=3,support_size=len(coefficients),l1_norm=2880,
              terms=[dict(word=[x+1 for x in w],coefficient=c)
                     for w,c in sorted(coefficients.items())])
    Path(path).write_text(json.dumps(data,indent=2)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--export",type=Path)
    parser.add_argument("--report",type=Path)
    parser.add_argument("--samples",type=int,default=12)
    args=parser.parse_args()
    result=verify(args.samples)
    if args.export:
        export(args.export,generate_coefficients_capelli())
    if args.report:
        args.report.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
