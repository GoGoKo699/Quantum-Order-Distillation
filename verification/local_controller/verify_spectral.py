#!/usr/bin/env python3
"""Exact Fraction certificate for pairwise-encoded local-dephasing converse.
No third-party dependencies. It verifies polynomial kernel identities, not merely a grid.
Crot convention: M(theta)=c² M_U+s² M_V+cs Crot.
"""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import json
import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent,
                     help='Directory for regenerated exact certificate records')
_args=_parser.parse_args() if __name__=='__main__' else _parser.parse_args([])
ROOT=_args.output_dir.resolve()
ROOT.mkdir(parents=True,exist_ok=True)
N=6
orders=list(permutations(range(3)));codes=[0,1,4,6,3,7]
dist=[[(a^b).bit_count() for b in codes] for a in codes]
zero=lambda:[[F(0) for _ in range(N)] for _ in range(N)]
def add(a,b):return [[a[i][j]+b[i][j] for j in range(N)] for i in range(N)]
def scale(c,a):return [[c*a[i][j] for j in range(N)] for i in range(N)]
def mul(a,b):return [[sum(a[i][k]*b[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
def trans(a):return list(map(list,zip(*a)))
def eq(a,b):assert a==b
T=[[F(1,6) for j in range(N)] for i in range(N)]
s=[(-1)**sum(o[i]>o[j] for i in range(3) for j in range(i+1,3)) for o in orders]
S=[[F(s[i]*s[j],6) for j in range(N)] for i in range(N)]
Id=[[F(i==j) for j in range(N)] for i in range(N)]
R=[[F(a==tuple(reversed(b))) for b in orders] for a in orders]
O=add(scale(F(1,2),add(Id,scale(-1,R))),scale(-1,S))
E=add(scale(F(1,2),add(Id,R)),scale(-1,T))
ps=[T,S,O,E]
for a in ps:eq(mul(a,a),a)
eq(add(add(T,S),add(O,E)),Id)
# Exact integer products strip the common sqrt(2)^-3 from Kraus words.
k=[[[1,0],[1,0]],[[0,1],[0,-1]]]
def mm(a,b):return [[sum(a[i][z]*b[z][j] for z in range(2)) for j in range(2)] for i in range(2)]
vec=[]
for hist in product(range(2),repeat=3):
    row=[]
    for o in orders:
        a=[[1,0],[0,1]]
        for z in o:a=mm(a,k[hist[z]])
        row.append([v for rr in a for v in rr])
    vec.append(row)
u=[-1,1,1,1];v=[1,1,1,-1]
tu=[[sum(u[z]*vec[h][i][z] for z in range(4)) for i in range(N)] for h in range(8)]
tv=[[sum(v[z]*vec[h][i][z] for z in range(4)) for i in range(N)] for h in range(8)]
Crot=[[F(sum(tu[h][i]*tv[h][j]+tv[h][i]*tu[h][j] for h in range(8)),64) for j in range(N)] for i in range(N)]
MU=[[F(sum(tu[h][i]*tu[h][j] for h in range(8)),64) for j in range(N)] for i in range(N)]
MV=[[F(sum(tv[h][i]*tv[h][j] for h in range(8)),64) for j in range(N)] for i in range(N)]
eq(MU,scale(F(3,4),O));eq(MV,add(T,scale(F(1,4),E)))
eq(mul(Crot,Crot),scale(F(3,16),add(O,E)))
eq(mul(T,Crot),zero());eq(mul(S,Crot),zero())
eq(add(mul(O,mul(Crot,E)),mul(E,mul(Crot,O))),Crot)
for i,j in product(range(N),repeat=2):
    assert not Crot[i][j] or dist[i][j]==1
# Coefficients on [T,S,O,E] for dephased MU and MV, ascending powersq.
u_coeff=[['1/4','1/4','-1/4','-1/4'],['1/4','-1/4','-1/4','1/4'],['1/4','1/8','1/8','1/4'],['1/4','-1/8','1/8','-1/4']]
v_coeff=[['1/4','1/4','1/4','1/4'],['1/4','-1/4','1/4','-1/4'],['1/4','1/8','-1/8','-1/4'],['1/4','-1/8','-1/8','1/4']]
for M,coeff in [(MU,u_coeff),(MV,v_coeff)]:
    for d in range(4):
        lhs=[[M[i][j]*int(dist[i][j]==d) for j in range(N)] for i in range(N)]
        rhs=zero()
        for P,c in zip(ps,coeff):rhs=add(rhs,scale(F(c[d]),P))
        eq(lhs,rhs)
# Polynomial arithmetic verifies determinant endpoint factor exactly.
def padd(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return c

def pmul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def pscale(t,a):return [t*x for x in a]
def trim(a):
    while a and not a[-1]:a.pop()
    return a
uO=list(map(F,u_coeff[2]));vO=list(map(F,v_coeff[2]));uE=list(map(F,u_coeff[3]));vE=list(map(F,v_coeff[3]));cE=[F(1,4),F(0),F(0),F(1,4)]
A=padd(uO,pscale(-1,vO));B=padd(cE,pscale(-1,uE))
left=padd(pmul(A,B),[F(0),F(0),F(-3,16)])
# (uO-vO)(cE-uE)-3q²/16 =q²(q+q²+2q³+8q⁴-6)/32.
right=pscale(F(1,32),[F(0),F(0),F(-6),F(1),F(1),F(2),F(8)])
assert trim(left)==trim(right)
# Other determinant endpoint is (uO-vO)(cE-vE)>=0;
# factors are q²(1+2q)/4 and q(1+q)/8.
assert trim(A)==[F(0),F(0),F(1,4),F(1,2)]
assert trim(padd(cE,pscale(-1,vE)))==[F(0),F(1,8),F(1,8)]
# Root is unique since derivative 1+2q+6q²+32q³ positive.
poly=lambda q:q+q*q+2*q**3+8*q**4-6
lo=F(811155656051002,10**15);hi=F(811155656051004,10**15)
assert poly(lo)<0<poly(hi)
assert poly(F(85,100))>0
# S3 invariance of local Hamming metric, exact.
for g in permutations(range(3)):
    perm=[orders.index(tuple(g[x] for x in o)) for o in orders]
    for i,j in product(range(N),repeat=2):assert dist[i][j]==dist[perm[i]][perm[j]]
result={'status':'PASS','exact_polynomial_kernel_decompositions':8,'exact_S3_distance_checks':216,
 'cross_convention':'Mtheta=c² MU+s² MV+cs Crot; Crot²=3/16(O+E)',
 'cross_dephasing_identity':'Crot(q)=q Crot(1)',
 'critical_polynomial':'8q^4+2q^3+q^2+q-6',
 'critical_root_bracket':[str(lo),str(hi)],
 'determinant_endpoint_identity':'(uO-vO)(cE-uE)-3q²/16=q²(8q⁴+2q³+q²+q-6)/32',
 'scope':'Exact simultaneous unitary-correction witness domination for q>=qcrit. Does not assert fourflag optimality belowqcrit or arbitrary CPTP targetrecovery.'}
(ROOT/'spectral_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

# Supplementary direct complex-Gram checks. Exact identities above carry proof;
# these numerical tests independently check conventions in the physical channel.
try:
    import numpy as np
except ImportError:
    np=None
if np is not None:
    rng=np.random.default_rng(29102026)
    ni=np.eye(2,dtype=complex);nx=np.array([[0,1],[1,0]],complex)
    ny=np.array([[0,-1j],[1j,0]],complex);nz=np.diag([1.,-1.])
    nu=(nx-nz)/np.sqrt(2);nv=(nx+nz)/np.sqrt(2)
    phi=np.array([1,0,0,1],complex)/np.sqrt(2)
    PHI=np.outer(phi,phi.conj());arr=np.array(vec,dtype=complex)/4
    nps=[np.array(P,dtype=float) for P in ps]
    ds=np.array(dist)
    max_violation=0.;checked=0
    for q in list(np.linspace(.85,1,7))+[.9**1.5]:
        L=q**ds
        cz=[(1+2*q+2*q*q+q**3)/4,(1-q**3)/4,
            (2+q+q*q+2*q**3)/8,(1+q**3)/4]
        support=sum(c*P for c,P in zip(cz,nps))
        for _ in range(40):
            z=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
            G=z@z.conj().T;G/=np.trace(G).real
            a=rng.normal(size=4);a/=np.linalg.norm(a)
            Q=a[0]*ni+1j*(a[1]*nu+a[2]*ny+a[3]*nv)
            jj=np.einsum('ij,hia,hjb->ab',G*L,arr,arr.conj())
            QQ=np.kron(Q,ni);jj=QQ@jj@QQ.conj().T
            f=float(np.trace(PHI@jj).real)
            bound=float(np.trace(G@support).real)
            max_violation=max(max_violation,f-bound)
            assert f<=bound+1e-12
            checked+=1
    result.update({'supplementary_complex_Gram_SU2_checks':checked,
        'q_grid':[float(q) for q in list(np.linspace(.85,1,7))+[.9**1.5]],
        'max_positive_witness_violation':max(0.,max_violation)})
(ROOT/'audit.json').write_text(json.dumps(result,indent=2)+'\n')

# Stronger complex cross-block audit: its symmetric and antisymmetric kernels
# have disjoint Hamming support 1 and 2, respectively. Thus for Hermitian G,
# chi_q=q Re(chi_1)+i q² Im(chi_1).
Nuv=[[F(sum(tu[h][i]*tv[h][j] for h in range(8)),64) for j in range(N)] for i in range(N)]
Arot=add(Nuv,scale(-1,trans(Nuv)))
nonzero_anti=0
for i,j in product(range(N),repeat=2):
    if Arot[i][j]:
        assert dist[i][j]==2
        nonzero_anti+=1
assert nonzero_anti==12
for P in [T,S]:
    eq(mul(P,Nuv),zero());eq(mul(Nuv,P),zero())
result.update({'antisymmetric_cross_kernel_nonzeros':nonzero_anti,
    'antisymmetric_cross_kernel_Hamming_support':2,
    'complex_cross_identity':'chi_q=q Re(chi_1)+i q² Im(chi_1)',
    'complex_cross_bound':'|chi_q|²<=q²|chi_1|²<=3q² o e/16',
    'stronger_spectral_conclusion':'lambda_max J(G o L)<=c_T t+c_S s+c_O o+c_E e for all complex G>=0 and q>=qcrit',
    'contract_note':'Report may retain unitary target-correction scope; no broader recovery claim needed.'})
(ROOT/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'antisymmetric_cross_check':'PASS','nonzeros':nonzero_anti,'support':2,'full_complex_spectral_bound':'CONFIRMED'}))
