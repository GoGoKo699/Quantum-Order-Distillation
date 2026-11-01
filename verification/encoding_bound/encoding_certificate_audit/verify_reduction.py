"""Independent exact audit of controller/correction/orbit reductions.

Standard library only. This does not replace the per-encoding PSD certificate.
"""
from fractions import Fraction as F
from itertools import permutations,product,combinations
from pathlib import Path
import json

import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
_args=_parser.parse_args()
HERE=_args.output_dir.resolve()
HERE.mkdir(parents=True,exist_ok=True)
orders=list(permutations(range(3)))
K=[[[1,0],[1,0]],[[0,1],[0,-1]]]
def mm(a,b):return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
vec=[]
for h in product(range(2),repeat=3):
    row=[]
    for o in orders:
        a=[[1,0],[0,1]]
        for j in o:a=mm(a,K[h[j]])
        row.append([v for rr in a for v in rr])
    vec.append(row)
basis=[[1,0,0,1],[-1,1,1,1],[0,-1,1,0],[1,1,1,-1]]
den=[32,64,32,64]
traces=[[[sum(b[a]*v[i][a] for a in range(4)) for i in range(6)] for v in vec] for b in basis]
M=[[[F(sum(t[h][i]*t[h][j] for h in range(8)),d) for j in range(6)] for i in range(6)] for t,d in zip(traces,den)]
H=[[F(sum(sum(v[i][a]*v[j][a] for a in range(4)) for v in vec),16) for j in range(6)] for i in range(6)]
N=[[F(sum(traces[1][h][i]*traces[3][h][j] for h in range(8)),64) for j in range(6)] for i in range(6)]
offdiag=0
for a,b in combinations(range(4),2):
    if (a,b)==(1,3):continue
    for i,j in product(range(6),repeat=2):
        assert sum(traces[a][h][i]*traces[b][h][j] for h in range(8))==0
        offdiag+=1
# Slot relabeling leaves every kernel, including the nonsymmetric UV kernel,
# exactly invariant. A code assignment permuted in the same way therefore has
# the same feasible support value, without a uniform-controller assumption.
invariance=0
for g in permutations(range(3)):
    p=[orders.index(tuple(g[x] for x in o)) for o in orders]
    for A in [H,*M,N]:
        for i,j in product(range(6),repeat=2):
            assert A[i][j]==A[p[i]][p[j]]
            invariance+=1
# The witness normalization is fixed by independent direct traces.
for i,j in product(range(6),repeat=2):
    assert H[i][j]==sum(A[i][j] for A in M)
# Exact local visibility field and outward rational bracket.
slo=F('3.1622776601683793319988935444327185337195551393')
shi=F('3.1622776601683793319988935444327185337195551394')
assert slo*slo<10<shi*shi
qlo,qhi=27*slo/100,27*shi/100
assert qlo*qlo<F(729,1000)<qhi*qhi
# Explicit numerical-objective constants are rational, not rounded inputs.
B=F(69,4);c=1-B/2
assert c==F(-61,8)
assert 1+B*F(1,8)==F(101,32)
result={'status':'PASS','zero_non_UV_Choi_kernel_entries':offdiag,
 'exact_S3_kernel_invariances':invariance,
 'trace_kernel_equals_sum_Bell_diagonal_kernels':True,
 'support_B':str(B),'support_c':str(c),'benchmark_cost':str(F(101,32)),
 'q_field':'27*sqrt(10)/100','q_squared':'729/1000',
 'sqrt10_bracket':[str(slo),str(shi)],'q_bracket':[str(qlo),str(qhi)],
 'scope':'Kernel/correction and physical-slot quotient audit; per-encoding PSD feasibility must be separately verified.'}
(HERE/'reduction_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
