"""Independent exact one/two-call support and attaining baseline certificate.

Standard library only. This derives the two-call Bell kernels directly from
integer Kraus words and proves the preparation-independent support witness.
"""
from fractions import Fraction as F
from itertools import product,combinations
from pathlib import Path
import json

import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
_args=_parser.parse_args()
HERE=_args.output_dir.resolve()
HERE.mkdir(parents=True,exist_ok=True)
KINT=[[[1,0],[1,0]],[[0,1],[0,-1]]]
def mm(a,b):return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
vec=[]
for hist in product(range(2),repeat=2):
    row=[]
    for order in [(0,1),(1,0)]:
        a=[[1,0],[0,1]]
        for k in order:a=mm(a,KINT[hist[k]])
        row.append([v for rr in a for v in rr])
    vec.append(row)
basis=[[1,0,0,1],[-1,1,1,1],[0,-1,1,0],[1,1,1,-1]]
traces=[[[sum(b[a]*v[i][a] for a in range(4)) for i in range(2)] for v in vec] for b in basis]
M=[]
for tr,den in zip(traces,[16,32,16,32]):
    M.append([[F(sum(tr[h][i]*tr[h][j] for h in range(4)),den) for j in range(2)] for i in range(2)])
H=[[F(sum(sum(v[i][a]*v[j][a] for a in range(4)) for v in vec),8) for j in range(2)] for i in range(2)]
assert H==[[F(1),F(1,2)],[F(1,2),F(1)]]
plus=[[F(1,4),F(1,4)],[F(1,4),F(1,4)]]
minus=[[F(1,4),F(-1,4)],[F(-1,4),F(1,4)]]
assert M==[plus,minus,plus,plus]
crosschecks=0
for a,b in combinations(range(4),2):
    for i,j in product(range(2),repeat=2):
        assert sum(traces[a][h][i]*traces[b][h][j] for h in range(4))==0
        crosschecks+=1
x=F(9,10);B=F(69,4);c=F(-61,8)
L=[[F(1),x],[x,F(1)]]
Hx=[[H[i][j]*L[i][j] for j in range(2)] for i in range(2)]
Ws=[[[L[i][j]*(B*A[i][j]+c*H[i][j]) for j in range(2)] for i in range(2)] for A in M]
WI,WU,_,_=Ws
K2=[[F(2),F(-2)],[F(-2),F(2)]]
def eig_pm(A):
    assert A[0][0]==A[1][1] and A[0][1]==A[1][0]
    return [A[0][0]+A[0][1],A[0][0]-A[0][1]]
assert eig_pm(WI)==[F(-229,80),F(-301,80)]
assert eig_pm(WU)==[F(-85,8),F(4)]
for W in Ws:
    slack=[[K2[i][j]-W[i][j] for j in range(2)] for i in range(2)]
    assert min(eig_pm(slack))>=0
assert eig_pm(K2)==[0,4]
# Native acceptance of minus has coefficient matrix Pminus/2.
Gminus=[[F(1,4),F(-1,4)],[F(-1,4),F(1,4)]]
trace=lambda A,G:sum(A[i][j]*G[j][i] for i,j in product(range(2),repeat=2))
p2=trace(Hx,Gminus)
f2=trace([[M[1][i][j]*L[i][j] for j in range(2)] for i in range(2)],Gminus)
assert (p2,f2)==(F(11,40),F(19,80))
assert p2+B*(f2-p2/2)==2
# One call followed by Hadamard gives equal I,Z unitary mixture.
# Its Choi largest eigenvalue is 1/2, so all unitary corrections satisfy f<=p/2.
raw1,raw2=F(21,61),F(40,61)
P=raw1+raw2*p2
Fbell=raw1/2+raw2*f2
calls=raw1+2*raw2
assert P==F(32,61) and Fbell/P==F(5,8)
assert calls/P==F(101,32)
result={'status':'PASS','exact_two_call_cross_kernel_zeros':crosschecks,
 'H_x':'[[1,x/2],[x/2,1]]','M_I_Y_V':'[[1,x],[x,1]]/4',
 'M_U':'[[1,-x],[-x,1]]/4','K2':'[[2,-2],[-2,2]]=4 Pminus',
 'K2_diagonal_ceiling':'2','WI_eigenvalues_plus_minus':list(map(str,eig_pm(WI))),
 'WU_eigenvalues_plus_minus':list(map(str,eig_pm(WU))),
 'two_minus_probability':str(p2),'two_minus_Bell_numerator':str(f2),
 'raw_trial_one_two_probabilities':[str(raw1),str(raw2)],'acceptance':str(P),
 'accepted_Bell_fidelity':str(Fbell/P),'half_diamond_error':'3/8',
 'expected_calls_per_acceptance':str(calls/P),
 'scope':'One/two-call support converse covers arbitrary preparation, ancillas and unitary target correction; matching protocol attains baseline.'}
(HERE/'baseline_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
