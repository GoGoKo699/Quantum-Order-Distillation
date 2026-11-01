#!/usr/bin/env python3
"""Exact, standard-library-only global encoding certificate verifier.

No SDP solver or numerical eigenvalue is trusted. A rational matrix is proved
positive by a rational Gram factor plus a strictly diagonally dominant residual.
An exact visibility interval covers the algebraic physical visibility.
"""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
from math import lcm
import json,hashlib,time

ROOT=Path(__file__).resolve().parent.parent
HERE=Path(__file__).resolve().parent
QHAT=F(170763,200000)
DELTA=F(1,10**7)
PERTURBATION=F(3,100000)
B=F(69,4);C=F(-61,8)

def zeros(n):return [[F(0) for _ in range(n)] for _ in range(n)]
def mm2(a,b):return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]

def exact_kernels():
    orders=list(permutations(range(3)))
    ks=[[[1,0],[1,0]],[[0,1],[0,-1]]]
    vec=[]
    for h in product(range(2),repeat=3):
        row=[]
        for order in orders:
            a=[[1,0],[0,1]]
            for k in order:a=mm2(a,ks[h[k]])
            row.append([x for rr in a for x in rr])
        vec.append(row)
    reps=[[1,0,0,1],[-1,1,1,1],[0,-1,1,0],[1,1,1,-1]]
    traces=[[[sum(rep[z]*vec[h][i][z] for z in range(4)) for i in range(6)] for h in range(8)] for rep in reps]
    ms=[]
    for tr,den in zip(traces,[32,64,32,64]):
        ms.append([[F(sum(tr[h][i]*tr[h][j] for h in range(8)),den) for j in range(6)] for i in range(6)])
    for a,b in product(range(4),repeat=2):
        if a==b or {a,b}=={1,3}:continue
        assert all(sum(traces[a][h][i]*traces[b][h][j] for h in range(8))==0 for i,j in product(range(6),repeat=2))
    hmat=[[F(sum(sum(v[i][a]*v[j][a] for a in range(4)) for v in vec),16) for j in range(6)] for i in range(6)]
    cross=[[F(sum(traces[1][h][i]*traces[3][h][j]+traces[3][h][i]*traces[1][h][j] for h in range(8)),64) for j in range(6)] for i in range(6)]
    mi,mu,my,mv=ms
    def comb(a,b,ma,mb):return [[a*ma[i][j]+b*mb[i][j] for j in range(6)] for i in range(6)]
    wi=comb(B,C,mi,hmat);wy=comb(B,C,my,hmat)
    w0=[[B*(mu[i][j]+mv[i][j])/2+C*hmat[i][j] for j in range(6)] for i in range(6)]
    wc=[[B*(mu[i][j]-mv[i][j])/2 for j in range(6)] for i in range(6)]
    ws=[[B*cross[i][j]/2 for j in range(6)] for i in range(6)]
    # Each varying entry of every real slack is one coefficient times q^d.
    # This covers the realification of A1=-Wc/2+i*Ws/2 as well.
    coeffs=[wi,wy,w0,[[x/2 for x in r] for r in wc],[[x/2 for x in r] for r in ws]]
    maximum=max(abs(x) for mat in coeffs for row in mat for x in row)
    assert maximum<F(4)
    # Independently check all 16 full Choi coordinate kernels under slot relabeling.
    invariances=0
    for g in permutations(range(3)):
        p=[orders.index(tuple(g[x] for x in o)) for o in orders]
        for a,b in product(range(4),repeat=2):
            mat=[[sum(v[i][a]*v[j][b] for v in vec) for j in range(6)] for i in range(6)]
            assert all(mat[i][j]==mat[p[i]][p[j]] for i,j in product(range(6),repeat=2))
            invariances+=36
    return (wi,wy,w0,wc,ws),{'full_Choi_kernel_entry_invariances':invariances,'max_variable_coefficient':str(maximum),'forbidden_ordered_cross_kernel_identities':10}

def matrices(cert,kernels):
    codes=cert['codes'];scale=cert['scale']
    assert len(codes)==6 and len(set(codes))==6 and all(isinstance(c,int) and 0<=c<8 for c in codes)
    k=[[F(x,scale) for x in row] for row in cert['K_integer']]
    xr=[[F(x,scale) for x in row] for row in cert['X_real_integer']]
    xi=[[F(x,scale) for x in row] for row in cert['X_imag_integer']]
    assert all(k[i][j]==k[j][i] and xr[i][j]==xr[j][i] and xi[i][j]==-xi[j][i] for i,j in product(range(6),repeat=2))
    assert all(k[i][i]<=F(27,10) for i in range(6))
    ell=[[QHAT**((a^b).bit_count()) for b in codes] for a in codes]
    wi,wy,w0,wc,ws=[[[ell[i][j]*mat[i][j] for j in range(6)] for i in range(6)] for mat in kernels]
    blocks={'K':k,'I':[[k[i][j]-wi[i][j] for j in range(6)] for i in range(6)],
            'Y':[[k[i][j]-wy[i][j] for j in range(6)] for i in range(6)]}
    re=zeros(12);im=zeros(12)
    for i,j in product(range(6),repeat=2):
        ar=-wc[i][j]/2;ai=ws[i][j]/2
        re[i][j]=xr[i][j];im[i][j]=xi[i][j]
        re[i][j+6]=ar;im[i][j+6]=ai
        re[j+6][i]=ar;im[j+6][i]=-ai
        re[i+6][j+6]=k[i][j]-w0[i][j]-xr[i][j]
        im[i+6][j+6]=-xi[i][j]
    realified=zeros(24)
    for i,j in product(range(12),repeat=2):
        realified[i][j]=re[i][j];realified[i][j+12]=-im[i][j]
        realified[i+12][j]=im[i][j];realified[i+12][j+12]=re[i][j]
    assert all(realified[i][j]==realified[j][i] for i,j in product(range(24),repeat=2))
    blocks['continuum']=realified
    return blocks

def verify_gram_residual(mat,factor):
    n=len(mat);s=factor['scale'];g=factor['integer']
    assert len(g)==n and all(len(row)==n for row in g)
    assert all(isinstance(x,int) for row in g for x in row)
    den=s*s
    for row in mat:
        for x in row:den=lcm(den,x.denominator)
    gram_factor=den//(s*s)
    residual=[[mat[i][j].numerator*(den//mat[i][j].denominator)-gram_factor*sum(g[i][z]*g[j][z] for z in range(n)) for j in range(n)] for i in range(n)]
    margins=[residual[i][i]-sum(abs(residual[i][j]) for j in range(n) if j!=i) for i in range(n)]
    minimum=min(margins)
    # Mhat=GG^T+R, R>=minimum/den I by exact diagonal dominance.
    # True visibility changes every eigenvalue by less than PERTURBATION.
    assert minimum*PERTURBATION.denominator>PERTURBATION.numerator*den
    return F(minimum,den)

def canonical(codes):
    orders=list(permutations(range(3)))
    images=[]
    for g in permutations(range(3)):
        p=[orders.index(tuple(g[x] for x in o)) for o in orders]
        z=[codes[j]^codes[p[0]] for j in p]
        for perm in permutations(range(3)):
            images.append(tuple(sum(((c>>perm[j])&1)<<j for j in range(3)) for c in z))
    return min(images)

def verify_coverage(candidates):
    covered={canonical(c['codes']) for c in candidates}
    assert len(candidates)==len(covered)==97
    histogram={rep:0 for rep in covered}
    for codes in permutations(range(8),6):
        rep=canonical(codes)
        assert rep in covered
        histogram[rep]+=1
    assert sum(histogram.values())==20160
    sizes={str(n):list(histogram.values()).count(n) for n in sorted(set(histogram.values()))}
    assert sizes=={'48':3,'144':49,'288':45}
    return {'all_injective_assignments':20160,'combined_classes':97,'orbit_sizes':sizes}

def verify_baseline():
    """Exact one/two-call support witnesses and the attaining fresh-trial mix."""
    k=[[[1,0],[1,0]],[[0,1],[0,-1]]]
    traces=[[],[],[],[]]
    reps=[[1,0,0,1],[-1,1,1,1],[0,-1,1,0],[1,1,1,-1]]
    words=[]
    for h in product(range(2),repeat=2):
        row=[]
        for order in [(0,1),(1,0)]:
            a=mm2(k[h[order[0]]],k[h[order[1]]])
            row.append([v for rr in a for v in rr])
        words.append(row)
    x=F(9,10);ell=[[F(1),x],[x,F(1)]]
    h2=[[F(sum(sum(v[i][a]*v[j][a] for a in range(4)) for v in words),8)*ell[i][j] for j in range(2)] for i in range(2)]
    assert h2==[[F(1),x/2],[x/2,F(1)]]
    for r,rep in enumerate(reps):
        traces[r]=[[sum(rep[a]*v[i][a] for a in range(4)) for i in range(2)] for v in words]
    for a,b in product(range(4),repeat=2):
        if a!=b:
            assert all(sum(traces[a][h][i]*traces[b][h][j] for h in range(4))==0 for i,j in product(range(2),repeat=2))
    for r,den in enumerate([16,32,16,32]):
        mat=[[F(sum(traces[r][h][i]*traces[r][h][j] for h in range(4)),den)*ell[i][j] for j in range(2)] for i in range(2)]
        sign=-1 if r==1 else 1
        assert mat==[[F(1,4),sign*x/4],[sign*x/4,F(1,4)]]
        # These commuting kernels have eigenvalues (1+x)/4 and (1-x)/4.
        assert (1+x)/4==F(19,40)
    pminus=[[F(1,2),F(-1,2)],[F(-1,2),F(1,2)]]
    pplus=[[F(1,2),F(1,2)],[F(1,2),F(1,2)]]
    upper=[[B*F(19,40)*int(i==j)+C*h2[i][j] for j in range(2)] for i in range(2)]
    assert upper==[[4*pminus[i][j]-F(229,80)*pplus[i][j] for j in range(2)] for i in range(2)]
    assert all(4*pminus[i][i]==2 for i in range(2))
    # One-call Choi is a PSD Gram matrix with J^2=J/2, trace one.
    vecs=[[1,0,1,0],[0,1,0,-1]]
    j1=[[F(sum(v[i]*v[j] for v in vecs),4) for j in range(4)] for i in range(4)]
    assert sum(j1[i][i] for i in range(4))==1
    assert all(sum(j1[i][z]*j1[z][j] for z in range(4))==j1[i][j]/2 for i,j in product(range(4),repeat=2))
    t1,t2=F(21,61),F(40,61);pm,fm=F(11,40),F(19,80)
    p=t1+t2*pm;f=t1/2+t2*fm;calls=t1+2*t2
    assert p==F(32,61) and f==F(20,61) and calls==F(101,61)
    assert f/p==F(5,8) and calls/p==F(101,32)
    return {'one_call_spectral_maximum':'1/2','two_call_diagonal_support_cap':'2',
            'raw_trial_probabilities':['21/61','40/61'],'acceptance':'32/61',
            'accepted_Bell_fidelity':'5/8','half_diamond_error':'3/8','cost':'101/32'}

def verify(candidates_path=None,factors_path=None,output_dir=None):
    started=time.time()
    candidates_path=Path(candidates_path or ROOT/'encoding_dual/certificate_candidates.json')
    factors_path=Path(factors_path or HERE/'rational_gram_factors.json')
    candidates=json.loads(candidates_path.read_text())['certificates']
    factors=json.loads(factors_path.read_text())['certificates']
    assert len(candidates)==len(factors)==97
    # Exact enclosing interval for the positive algebraic q=sqrt(729/1000).
    assert 0<QHAT-DELTA<QHAT<1
    assert (QHAT-DELTA)**2<F(729,1000)<QHAT**2
    # |q^d-qhat^d| <= 3 DELTA, coefficient<4, dimension<=24.
    assert 24*4*3*DELTA<PERTURBATION
    kernels,kernel_status=exact_kernels()
    minimum=None
    for index,(candidate,factor) in enumerate(zip(candidates,factors)):
        assert candidate['codes']==factor['codes']
        mats=matrices(candidate,kernels)
        assert set(mats)==set(factor['blocks'])
        for name,mat in mats.items():
            margin=verify_gram_residual(mat,factor['blocks'][name])
            minimum=margin if minimum is None else min(minimum,margin)
    coverage=verify_coverage(candidates)
    baseline=verify_baseline()
    result={'status':'PASS','certificate_count':97,'exact_positive_matrix_checks':388,
            'common_support_bound':'27/10','required_advantage_threshold':'3',
            'q_interval':[str(QHAT-DELTA),str(QHAT)],'q_squared':'729/1000',
            'uniform_operator_perturbation_bound':str(PERTURBATION),
            'minimum_exact_DD_margin':str(minimum),'minimum_exact_DD_margin_decimal':float(minimum),
            'kernel_checks':kernel_status,'coverage':coverage,'baseline':baseline,
            'candidate_SHA256':hashlib.sha256(candidates_path.read_bytes()).hexdigest(),
            'factor_SHA256':hashlib.sha256(factors_path.read_bytes()).hexdigest(),
            'runtime_seconds':round(time.time()-started,3),
            'conclusion':'Exact C123(3/8)=101/32 for injective six-label computational-basis assignments; final target unitary corrections; declared local phase noise.'}
    output=Path(output_dir) if output_dir is not None else HERE
    output.mkdir(parents=True,exist_ok=True)
    (output/'exact_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,help='Directory for the regenerated exact result')
    args=parser.parse_args()
    verify(output_dir=args.output_dir)
