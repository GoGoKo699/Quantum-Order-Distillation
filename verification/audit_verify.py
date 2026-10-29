#!/usr/bin/env python3
"""Independent exact verifier for the tilt certificates (Python stdlib + NumPy).

Default verifies the stored exact identities, an independent rational Kraus
point, all CT/RT grid records, and the full orthogonal-family coefficients. --rebuild also regenerates every four-call
polynomial coefficient and bounds the original int64 intermediates.
Expected adjacent files: tilt_exact_polynomials.npz and
joint_tilt_certificate_results.json, and family8_polynomials.npz.
Use --root DIR for a certificate bundle
or a source-tree root. No writes are performed.
This does not import the research scripts or write/replace their certificates.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import comb, gcd, atan
from pathlib import Path
import json
import argparse
import numpy as np


def certificate(root, name, original_subdir):
    direct = root/name
    return direct if direct.exists() else root/original_subdir/name


ORDERS = [tuple(int(c)-1 for c in w) for w in
          ('1234','1243','2134','2143','3412','3421','4312','4321')]
CODES = (0,1,2,3,4,6,5,7)
WALSH = np.array([[(-1)**((x&k).bit_count()) for k in range(8)]
                  for x in CODES], dtype=object)
GN = np.array([[4**((x^y).bit_count()) * 5**(3-(x^y).bit_count())
                for y in CODES] for x in CODES], dtype=object)
BASES = [np.array(a, dtype=object) for a in
         ([[1,0],[0,1]], [[0,1],[-1,0]], [[0,1],[1,0]], [[1,0],[0,-1]])]


def polynomial_products():
    """Rebuild with Python integers; independently bound every int64 sum."""
    factors = [
        [np.array([[1,0],[1,0]], dtype=object),
         np.array([[-1,0],[1,0]], dtype=object)],
        [np.array([[0,-1],[0,1]], dtype=object),
         np.array([[0,-1],[0,-1]], dtype=object)],
    ]
    for i,j in product(range(2), repeat=2):
        a = np.zeros((2,2), dtype=object); a[i,j] = 1
        factors.append([a])
    gg = np.zeros((9,4,4,8,8), dtype=object)
    bb = np.zeros((9,8,8), dtype=object)
    absolute_g = np.zeros_like(gg); absolute_b = np.zeros_like(bb)
    largest_v = largest_trace = 0
    for labs in product(range(6), repeat=4):
        k = sum(j < 2 for j in labs)
        vv = np.zeros((8,k+1,2,2), dtype=object)
        for o,order in enumerate(ORDERS):
            v = [np.eye(2, dtype=object)]
            for j in order:
                w = [np.zeros((2,2), dtype=object)
                     for _ in range(len(v)+len(factors[labs[j]])-1)]
                for p,a in enumerate(v):
                    for q,b in enumerate(factors[labs[j]]):
                        w[p+q] += a @ b
                v = w
            vv[o] = v
        largest_v = max(largest_v, max(abs(int(x)) for x in vv.flat))
        ts = np.array([[[sum((u @ vv[o,p]).diagonal()) for o in range(8)]
                        for p in range(k+1)] for u in BASES], dtype=object)
        largest_trace = max(largest_trace, max(abs(int(x)) for x in ts.flat))
        for p,q in product(range(k+1), repeat=2):
            term_g = np.multiply.outer(ts[:,p], ts[:,q]).transpose(0,2,1,3)
            term_b = vv[:,p].reshape(8,4) @ vv[:,q].reshape(8,4).T
            for r in range(5-k):
                scale = 99**k * comb(4-k,r)
                d = p+q+2*r
                gg[d] += scale*term_g; bb[d] += scale*term_b
                absolute_g[d] += scale*abs(term_g)
                absolute_b[d] += scale*abs(term_b)
    # Original einsum multiplies two +/-1 Walsh matrices. Bounding by the
    # absolute term sum * max(GN) * 8^2 also bounds all its partial sums.
    bound_g = int(np.max(absolute_g))*125*64
    bound_b = int(np.max(absolute_b))*125*64
    assert max(bound_g,bound_b) < 2**63
    gg *= GN; bb *= GN
    gx = np.empty_like(gg); bx = np.empty_like(bb)
    for p in range(9):
        bx[p] = WALSH.T @ bb[p] @ WALSH
        for i,j in product(range(4), repeat=2):
            gx[p,i,j] = WALSH.T @ gg[p,i,j] @ WALSH
    print('Python-integer coefficient reconstruction: PASS')
    print('int64 absolute intermediate bounds G/B:', bound_g, bound_b,
          '; product/trace coefficient maxima:', largest_v, largest_trace)
    return gx,bx


def convolution(a,b):
    c = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j] += int(x)*int(y)
    return c


def bernstein(poly, halfwidth):
    n = len(poly)-1; a = [F(0)]*(n+1)
    for j,c in enumerate(poly):
        for k in range(j+1):
            a[k] += c*comb(j,k)*halfwidth**j * 2**k * (-1)**(j-k)
    return [sum(a[j]*F(comb(k,j),comb(n,j)) for j in range(k+1))
            for k in range(n+1)]


def audit_blocks(gx,bx,stored):
    def diagonal(a):
        return all(a[p,i,j] == 0 for p,i,j in product(range(9),range(8),range(8))
                   if i != j)
    assert diagonal(bx)
    # Real Kraus matrices alone do not remove imaginary off-diagonal control
    # objectives for a complex controller effect. These exact cross-sector
    # symmetries are needed for the all-SU(2) convex-circle reduction.
    for i,j in product((0,1),(2,3)):
        assert np.array_equal(gx[:,i,j],gx[:,i,j].transpose(0,2,1))
        assert not np.any(np.diagonal(gx[:,i,j],axis1=1,axis2=2))
    first = (gx[:,0,0]+gx[:,1,1], gx[:,0,0]-gx[:,1,1],
             gx[:,0,1]+gx[:,1,0])
    c0 = gx[:,2,2]+gx[:,3,3]
    cc = gx[:,2,2]-gx[:,3,3]; cs = gx[:,2,3]+gx[:,3,2]
    assert all(diagonal(a) for a in (*first,c0))
    allowed = {(1,5),(5,1),(2,6),(6,2)}
    assert all(a[p,i,j] == 0 for a in (cc,cs)
               for p,i,j in product(range(9),range(8),range(8))
               if i != j and (i,j) not in allowed)
    h1,h2,k = (cc[:,1,1],cs[:,1,1]), (cc[:,5,5],cs[:,5,5]), (cc[:,1,5],cs[:,1,5])
    def dot(a,b):
        return np.array(convolution(a[0],b[0]), dtype=object) + convolution(a[1],b[1])
    assert np.array_equal(convolution(h1[0],h2[1]),convolution(h1[1],h2[0]))
    assert not np.any(dot(h1,k)); assert not np.any(dot(h2,k))
    for a in (cc,cs):
        assert np.array_equal(a[:,1,1],a[:,2,2])
        assert np.array_equal(a[:,5,5],a[:,6,6])
        assert np.array_equal(a[:,1,5],-a[:,2,6])
    margin = -dot(h1,h2)-dot(k,k)
    common = reduce(gcd, (abs(int(v)) for v in margin))
    primitive = [int(v)//common for v in margin]
    while primitive[-1] == 0: primitive.pop()
    assert primitive == list(stored['margin'])
    assert common == stored['margin_gcd'].item()
    b = bernstein(primitive,F(1,5))
    assert min(b) > 0
    print('Exact block identities and degree-%s Bernstein positivity: PASS' % (len(b)-1))
    print('minimum Bernstein coefficient:', min(b), '=', float(min(b)))


def direct_rational_point(gx,bx,t=F(1,10)):
    """No polynomial multiplication here: independently contract raw Kraus words."""
    # t=1/10 gives normalized pure Kraus column (9,11)/sqrt(202).
    m,n = 9,11; den=m*m+n*n
    ks = [np.array([[m,0],[n,0]],dtype=object),
          np.array([[0,-n],[0,m]],dtype=object)]
    for a,b in product(range(2),repeat=2):
        v=np.zeros((2,2),dtype=object); v[a,b]=1; ks.append(v)
    rawg=np.zeros((5,4,4,8,8),dtype=object)
    rawb=np.zeros((5,8,8),dtype=object)
    for labs in product(range(6),repeat=4):
        k=sum(j<2 for j in labs); vv=[]
        for order in ORDERS:
            a=np.eye(2,dtype=object)
            for j in order: a=a@ks[labs[j]]
            vv.append(a)
        vv=np.array(vv,dtype=object)
        tr=np.array([[sum((u@a).diagonal()) for a in vv] for u in BASES],dtype=object)
        rawg[k] += np.multiply.outer(tr,tr).transpose(0,2,1,3)
        rawb[k] += vv.reshape(8,4)@vv.reshape(8,4).T
    directg=np.zeros((4,4,8,8),dtype=object); directb=np.zeros((8,8),dtype=object)
    for k in range(5):
        weight=F(99,100)**k * F(1,200)**(4-k) / den**k
        directg += weight*rawg[k]*GN/F(4*125)
        directb += weight*rawb[k]*GN/F(2*125)
    directb=WALSH.T@directb@WALSH/F(8)
    for i,j in product(range(4),repeat=2):
        directg[i,j]=WALSH.T@directg[i,j]@WALSH/F(8)
    d=64*100**4*125*8*(1+t*t)**4
    pg=sum(gx[p]*t**p for p in range(9))/d
    pb=sum(bx[p]*t**p for p in range(9))/(d/2)
    assert np.array_equal(pg,directg); assert np.array_equal(pb,directb)
    assert sum(directb.diagonal()) == 8  # each normalized Walsh mode B; sum = r
    print('Independent t=1/10 Kraus contractions, G and B normalization: PASS')



# Independent CT/RT products use ordinary Python integer tuples, separately
# from the NumPy polynomial construction above.
I = (1, 0, 0, 1)
NOISE = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
ORD3 = [(0, 1, 2), (0, 2, 1), (1, 0, 2)]
ORD4 = ORDERS
SIGNS = [1, -1, -1, 1, 1, -1, -1, 1]
PAIRS = [(0, 1), (0, 2), (1, 2)]
DEPH = [[4**((a^b).bit_count()) * 5**(3-(a^b).bit_count())
         for b in CODES] for a in CODES]
CONTRACTION = [[SIGNS[i]*SIGNS[j]*DEPH[i][j] for j in range(8)] for i in range(8)]


def mm(A, B):
    a, b, c, d = A
    e, f, g, h = B
    return (a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h)


def ct_state(m, n):
    """CT output on |0>, balanced 3-order control, exactly normalized."""
    K = [(m, 0, n, 0), (0, -n, 0, m)] + NOISE
    raw = [[[0]*6 for _ in range(6)] for _ in range(4)]
    for labs in product(range(6), repeat=3):
        k = sum(x < 2 for x in labs)
        v = []
        for o in ORD3:
            a = mm(K[labs[o[2]]], mm(K[labs[o[1]]], K[labs[o[0]]]))
            v.extend((a[0], a[2]))
        for i in range(6):
            for j in range(6):
                raw[k][i][j] += v[i]*v[j]
    d = m*m+n*n
    return [[sum(F(raw[k][i][j]*99**k,
                   100**3 * 3 * d**k * 2**(3-k)) for k in range(4))
             for j in range(6)] for i in range(6)]


def ct_shifted_pt(r, pair):
    # Partial transpose of the two-dimensional target, then physical dephasing.
    return [[r[2*(i//2)+j%2][2*(j//2)+i%2]
             * F(4, 5)**(0 if i//2 == j//2 else
                          (2 if set((i//2, j//2)) == set(pair) else 1))
             - F(int(i == j), 100)
             for j in range(6)] for i in range(6)]


def ldl(A):
    L = [[F(int(i == j)) for j in range(6)] for i in range(6)]
    D = []
    for i in range(6):
        di = A[i][i] - sum(L[i][k]**2 * D[k] for k in range(i))
        assert di != 0
        D.append(di)
        for j in range(i+1, 6):
            L[j][i] = (A[j][i] - sum(L[j][k]*L[i][k]*D[k]
                                    for k in range(i))) / di
    # Explicit reconstruction, beyond positivity of the pivots.
    for i in range(6):
        for j in range(6):
            assert sum(L[i][k]*D[k]*L[j][k] for k in range(6)) == A[i][j]
    return D


def rt_values(m, n):
    """RT Bell-input p, unnormalized Bell overlap f, witness f-p/2.

    16000=64*125*2: balanced input and rank-one X outcome, dephasing,
    and normalized Choi trace.  The extra factor2 for f is Bell overlap.
    """
    K = [(m, 0, n, 0), (0, -n, 0, m)] + NOISE
    rawp, rawf = [0]*5, [0]*5
    for labs in product(range(6), repeat=4):
        k = sum(x < 2 for x in labs)
        vs = []
        for o in ORD4:
            a = I
            for z in o:
                a = mm(a, K[labs[z]])
            vs.append(a)
        tr = [v[0]+v[3] for v in vs]
        for i in range(8):
            for j in range(8):
                rawp[k] += CONTRACTION[i][j] * sum(
                    vs[i][z]*vs[j][z] for z in range(4))
                rawf[k] += CONTRACTION[i][j] * tr[i]*tr[j]
    d = m*m+n*n
    p = sum(F(rawp[k]*99**k, 100**4*d**k*2**(4-k)*16000)
            for k in range(5))
    f = sum(F(rawf[k]*99**k, 100**4*d**k*2**(4-k)*32000)
            for k in range(5))
    return p, f, f-p/2


def audit_ct_rt(root=None):
    root = Path(root) if root is not None else Path(__file__).resolve().parent
    records = json.loads(certificate(root, 'joint_tilt_certificate_results.json',
                                    'entanglement_bridge_analysis').read_text())
    ct = {(r['k'], tuple(r['long_pair'])): r for r in records['CT']}
    rt = {r['k']: r for r in records['RT']}
    assert len(ct) == len(records['CT']) == 249
    assert len(rt) == len(records['RT']) == 83
    assert set(ct) == {(k, p) for k in range(-41, 42) for p in PAIRS}
    assert set(rt) == set(range(-41, 42))
    minpivot = None
    selected = {}
    for k in range(-41, 42):
        r = ct_state(400-k, 400+k)
        assert sum(r[i][i] for i in range(6)) == 1
        assert r == [list(c) for c in zip(*r)]
        for pair in PAIRS:
            row = ct[k, pair]
            assert F(row['shift']) == F(1, 100)
            D = ldl(ct_shifted_pt(r, pair))
            assert D == [F(x) for x in row['ldl_pivots']]
            assert all(d > 0 for d in D)
            cand = (min(D), k, pair, D.index(min(D)))
            if minpivot is None or cand[0] < minpivot[0]:
                minpivot = cand
        p, f, w = rt_values(400-k, 400+k)
        assert (p, f, w) == tuple(F(rt[k][x]) for x in ('p', 'f', 'witness'))
        assert w > F(1, 125)
        assert F(1, 25) < p < F(9, 100)
        if k in (-41, 0, 41):
            selected[k] = {'p': float(p), 'F': float(f/p), 'w': float(w)}
    # These floating evaluations display coverage only.  Rigorous coverage:
    # adjacent parameter gaps <=2*(1/400)=1/200 by atan's 1-Lipschitz bound;
    # atan(41/400)>=41/400-(41/400)^3/3 >1/10 proves endpoint coverage.
    x = F(41, 400)
    assert 2*(x-x**3/3) > F(1, 5)
    # Claimed continuum constants from conservative nearest distance<=1/400.
    assert F(1,100)-3*F(1,400) == F(1,400)
    assert F(1,125)-2*F(1,400) == F(3,1000)
    assert F(1,25)-2*F(1,400) == F(7,200)
    assert F(9,100)+2*F(1,400) == F(19,200)
    assert F(1,2)+F(3,1000)/F(19,200) == F(101,190)
    return {'CT_exact_records_verified': len(ct),
            'RT_exact_records_verified': len(rt),
            'min_CT_LDL_pivot': float(minpivot[0]),
            'min_CT_LDL_location': [minpivot[1], list(minpivot[2]), minpivot[3]],
            'selected_RT': selected,
            'grid_endpoint': 2*atan(41/400),
            'grid_max_gap': 2*atan(1/400),
            'continuum_CT_strict_eigenvalue_lower_bound': '1/400',
            'continuum_RT_strict_witness_lower_bound': '3/1000',
            'continuum_RT_p_strict_bounds': ['7/200', '19/200'],
            'continuum_RT_F_strict_lower_bound': '101/190'}





WORDS = ('1234','1243','2134','2143','3412','3421','4312','4321')
ORDERS = [tuple(int(c)-1 for c in w) for w in WORDS]
CODES = (0,1,2,3,4,6,5,7)
I = (1,0,0,1)
L = [(1,0,1,0), (0,1,0,-1),
     (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)]
U = [I, (0,1,-1,0), (0,1,1,0), (1,0,0,-1)]  # I, J=iY, X, Z


def mm(a,b):
    a0,a1,a2,a3=a; b0,b1,b2,b3=b
    return (a0*b0+a1*b2, a0*b1+a1*b3,
            a2*b0+a3*b2, a2*b1+a3*b3)


def trace_product(a,b):
    return a[0]*b[0]+a[1]*b[2]+a[2]*b[1]+a[3]*b[3]


def coefficients():
    Bbern=np.zeros((5,8,8),dtype=object)
    Gbern=np.zeros((5,4,4,8,8),dtype=object)
    for labs in product(range(6),repeat=4):
        k=sum(a<2 for a in labs)
        vs=[]
        for o in ORDERS:
            a=I
            for z in o:
                a=mm(a,L[labs[z]])
            vs.append(a)
        ts=[[trace_product(u,v) for v in vs] for u in U]
        for i in range(8):
            for j in range(8):
                Bbern[k,i,j]+=sum(vs[i][z]*vs[j][z] for z in range(4))
                for a in range(4):
                    for b in range(4):
                        Gbern[k,a,b,i,j]+=ts[a][i]*ts[b][j]
    # The common physical product denominator is 2^4=16.  Normalized
    # Choi trace contributes 1/2 (B/32), Bell overlap 1/4 (A/64).
    Bp=np.zeros_like(Bbern); Gp=np.zeros_like(Gbern)
    for k in range(5):
        for j in range(5-k):
            c=(-1)**j*comb(4-k,j)
            Bp[k+j]+=c*Bbern[k]
            Gp[k+j]+=c*Gbern[k]
    # Endpoint trace numerators for (X+Z)/sqrt(2), (X-Z)/sqrt(2).
    plus=Gp[:,2,2]+Gp[:,3,3]+Gp[:,2,3]+Gp[:,3,2]
    minus=Gp[:,2,2]+Gp[:,3,3]-Gp[:,2,3]-Gp[:,3,2]
    assert all(int(x)%2==0 for x in plus.flat)
    assert all(int(x)%2==0 for x in minus.flat)
    A=np.stack((Gp[:,0,0],Gp[:,1,1],plus//2,minus//2))
    return Bp,Gp,A


def bernstein_unit_interval(p):
    n=len(p)-1
    return [sum(F(int(p[j]))*F(comb(k,j),comb(n,j))
                for j in range(k+1)) for k in range(n+1)]


def audit_orthogonal(root=None):
    root=Path(root) if root is not None else Path(__file__).resolve().parents[1]
    B,G,A=coefficients()
    source=root/'family8_polynomials.npz'
    if not source.exists():
        source=root/'fourcall_optimum_analysis'/'family8_polynomials.npz'
    assert source.exists(), source
    with np.load(source) as z:
        for name,ours in [('Bp',B),('Gp',G),('Aends',A)]:
            assert np.array_equal(ours,z[name]), name
        assert np.array_equal(z['orders'],ORDERS)
    H=np.array([[(-1)**((code&j).bit_count()) for j in range(8)]
                for code in CODES],dtype=object)
    assert np.array_equal(H.T@H,8*np.eye(8,dtype=object))
    dist=np.array([[(x^y).bit_count() for y in CODES]for x in CODES])
    masks=np.array([dist==h for h in range(4)],dtype=object)
    # Each independent (s^k gamma^h) coefficient, not a numerical sample.
    diagonal_checked=0
    for M in list(B)+list(A.reshape(-1,8,8)):
        for mask in masks:
            W=H.T@(M*mask)@H
            assert np.array_equal(W,np.diag(np.diag(W)))
            diagonal_checked+=1
    # Parallel-circle sine and perpendicular-circle cosine vanish identically.
    assert np.all(G[:,0,1]+G[:,1,0]==0)
    assert np.all(G[:,2,2]-G[:,3,3]==0)
    # Cross-sector symmetry eliminates imaginary interference for every SU(2) U.
    # Vanishing Walsh diagonal cross-sector coefficients gives scalar success
    # effects for each real physical Walsh branch (including local dephasing).
    for a in (0,1):
        for b in (2,3):
            for M in G[:,a,b]:
                assert np.array_equal(M,M.T)
                for mask in masks:
                    assert np.all(np.diag(H.T@(M*mask)@H)==0)
    local=root/'family8_local_polynomials.npz'
    if not local.exists():
        local=root/'fourcall_optimum_analysis'/'family8_local_polynomials.npz'
    if local.exists():
        with np.load(local) as z:
            assert np.array_equal(z['Bl'],B[:,None]*masks[None])
            assert np.array_equal(z['Al'],A[:,:,None]*masks[None,None])
            assert np.array_equal(z['dist'],dist)
            assert np.array_equal(z['codes'],CODES)
    b=np.zeros((8,5),dtype=object); a=np.zeros((8,4,5),dtype=object)
    for j in range(8):
        v=H[:,j]
        for k in range(5):
            assert (v@B[k]@v)%8==0
            b[j,k]=(v@B[k]@v)//8
            for end in range(4):
                assert (v@A[end,k]@v)%8==0
                a[j,end,k]=(v@A[end,k]@v)//8
    # Check all six distinct classes in orthogonal_section.tex.  The arrays
    # below are numerator polynomials in ascending powers of s: b/32,a/64.
    expected_b=[(76,0,32,0,4),(36,0,-16,0,-4),(36,0,-16,0,-4),
                (12,0,0,0,4),(24,0,8,0,0),(24,0,8,0,0),
                (24,0,8,0,0),(24,0,-24,0,0)]
    assert np.array_equal(b,np.array(expected_b,dtype=object))
    candidates=[(0,),(0,3),(0,3),(0,),(3,),(1,),(1,),(2,)]
    expected_best=[[(56,0,32,0,8)],[(24,0,-16,0,-8),(16,16,0,0,0)],
                   [(24,0,-16,0,-8),(16,16,0,0,0)],[(24,0,0,0,8)],
                   [(16,16,16,16,0)],[(16,0,16,0,0)],
                   [(16,0,16,0,0)],[(16,16,-16,-16,0)]]
    for j in range(8):
        assert [tuple(a[j,c]) for c in candidates[j]]==expected_best[j]
        for end in range(4):
            # Every other endpoint is below at least one claimed candidate
            # on the full interval s in [0,1], by exact Bernstein coefficients.
            assert any(all(c>=0 for c in bernstein_unit_interval(a[j,cand]-a[j,end]))
                       for cand in candidates[j])
    # N_s^4=Delta: B's branch diagonal is 1; all endpoint A diagonals 1/4.
    for i in range(8):
        assert tuple(B[:,i,i])==(32,0,0,0,0)
        for end in range(4):
            assert tuple(A[end,:,i,i])==(16,0,0,0,0)
    # Local physical example, calculated from rebuilt independent coefficients.
    s=F(99,100); gamma=F(4,5)
    def evaluate(M):
        return sum(M[k]*s**k for k in range(5)) * np.array(
            [[gamma**int(dist[i,j])for j in range(8)]for i in range(8)],dtype=object)
    Be=evaluate(B)/32
    Ae=[evaluate(A[end])/64 for end in range(4)]
    rows=[]
    for j in range(8):
        v=H[:,j]
        bj=(v@Be@v)/8
        rewards=[(v@M@v)/8 for M in Ae]
        aj=max(rewards)
        rows.append({'mode':j,'capacity':str(bj/8),'best_F':str(aj/bj),
                     'endpoint':rewards.index(aj)})
    rows.sort(key=lambda r:F(r['best_F']),reverse=True)
    assert F(rows[0]['best_F'])==F(1566163803,2065458803)
    assert F(rows[0]['capacity'])==F(6196376409,10**11)
    det=sum(F(r['capacity'])*F(r['best_F']) for r in rows)
    assert det==F(23664816609,50000000000)
    return {'independent_Kraus_tuples':6**4,
            'coefficient_matrices_Walsh_diagonal':diagonal_checked,
            's_polynomial_degree':4,'gamma_polynomial_degree':3,
            'source_arrays_exact_match':(['Bp','Gp','Aends','Bl','Al'] if local.exists()
                                         else ['Bp','Gp','Aends']),
            'zero_harmonics':True,'cross_sector_symmetry':True,
            'Walsh_cross_sector_diagonals_zero':True,
            'six_class_table_and_dominance_on_unit_interval':True,
            'B_denominator':32,'A_denominator':64,
            'local_example':rows,'local_Fdet':str(det)}




def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent,
                        help='directory containing the certificates (or source root)')
    parser.add_argument('--rebuild',action='store_true',
                        help='also regenerate all four-call polynomial coefficients')
    args=parser.parse_args()
    root=args.root.resolve()
    stored=np.load(certificate(root,'tilt_exact_polynomials.npz','tilt_analysis'),
                   allow_pickle=True)
    assert stored['Gx'].shape == (9,4,4,8,8)
    assert stored['Bx'].shape == (9,8,8)
    if args.rebuild:
        gx,bx=polynomial_products()
        assert np.array_equal(gx,stored['Gx']); assert np.array_equal(bx,stored['Bx'])
        print('All stored polynomial coefficients: PASS')
    else:
        gx=stored['Gx'].astype(object); bx=stored['Bx'].astype(object)
    assert np.array_equal(WALSH,stored['H'])
    print('Walsh matrix: PASS')
    audit_blocks(gx,bx,stored)
    direct_rational_point(gx,bx)
    # CT/RT section follows; it also uses a separate raw-Kraus construction.
    print(json.dumps(audit_ct_rt(root), indent=2))
    print(json.dumps(audit_orthogonal(root), indent=2))


if __name__ == '__main__':
    main()
