from fractions import Fraction as F
import itertools as it,json
from pathlib import Path
import argparse
_parser=argparse.ArgumentParser(description='Exact constructive alternative-encoding certificate')
_parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
_args=_parser.parse_args()
_output_dir=_args.output_dir.resolve()
_output_dir.mkdir(parents=True,exist_ok=True)
ORDERS=list(it.permutations(range(3)));CODES=[0,3,5,2,4,1]
t=F(49,50);q=t**3;qdist=[[q**((a^b).bit_count())for b in CODES]for a in CODES]
par=[(-1)**sum(o[i]>o[j]for i in range(3)for j in range(i+1,3))for o in ORDERS]
T=[[F(1,6)for _ in range(6)]for _ in range(6)]
S=[[F(par[i]*par[j],6)for j in range(6)]for i in range(6)]
R=[[int(ORDERS[i]==ORDERS[j][::-1])for j in range(6)]for i in range(6)]
O=[[F(int(i==j)-R[i][j],2)-S[i][j]for j in range(6)]for i in range(6)]
E=[[F(int(i==j)+R[i][j],2)-T[i][j]for j in range(6)]for i in range(6)]
M=[[[1,0],[1,0]],[[0,1],[0,-1]]]
def mm(a,b):return[[sum(a[i][k]*b[k][j]for k in range(2))for j in range(2)]for i in range(2)]
vecs=[]
for hist in it.product(range(2),repeat=3):
 row=[]
 for order in ORDERS:
  a=[[1,0],[0,1]]
  for slot in order:a=mm(a,M[hist[slot]])
  row.append([z for rr in a for z in rr])
 vecs.append(row)
Q=[[1,0,0,1],[-1,1,1,1],[0,-1,1,0],[1,1,1,-1]]
weights=[];ch=[]
for name,P in zip(['T','S','O','E'],[T,S,O,E]):
 J=[[sum(P[i][j]*qdist[i][j]*sum(v[i][a]*v[j][b]for v in vecs)for i,j in it.product(range(6),repeat=2))/96 for b in range(4)]for a in range(4)]
 w=[sum(z[a]*J[a][b]*z[b]for a,b in it.product(range(4),repeat=2))/den for z,den in zip(Q,[2,4,2,4])]
 # All offdiag between realBellIUV vanish. Ycross vanishes fromrealJ + oldkernelidentity.
 for a,b in it.combinations(range(4),2):assert sum(Q[a][i]*J[i][j]*Q[b][j]for i,j in it.product(range(4),repeat=2))==0
 assert sum(w)==sum(J[i][i]for i in range(4))
 for a,b in it.product(range(2),repeat=2):assert J[a][b]+J[a+2][b+2]==sum(w)*F(int(a==b),2)
 assert all(z>=0 for z in w)
 p=sum(w);f=max(w);weights.append((p,f));ch.append(J)
 print(name,'weights',list(map(str,w)),'p',p,'f',f,'error',1-f/p,'correction', ['I','U','Y','V'][w.index(f)])
# allsubsetpolicies +baseline vertices
verts=[(F(1,2),F(1),F(1),'1')]
x=t**2;p2=(2-x)/4;f2=(1+x)/8
verts.append((1-f2/p2,2/p2,p2,'2minus'))
verts.append((1-2*f2,F(2),F(1),'2both'))
for mask in range(1,16):
 ps=sum(weights[k][0]for k in range(4)if mask>>k&1);fs=sum(weights[k][1]for k in range(4)if mask>>k&1)
 verts.append((1-fs/ps,3/ps,ps,'3'+''.join(['T','S','O','E'][k]for k in range(4)if mask>>k&1)))
eps=F(3,8);best=None
for a,b in it.combinations(verts,2):
 if a[0]==b[0]:continue
 alpha=(eps-b[0])/(a[0]-b[0])
 if not 0<=alpha<=1:continue
 cost=alpha*a[1]+(1-alpha)*b[1]
 if best is None or cost<best[0]:best=(cost,alpha,a,b)
co,alpha,a,b=best
raw_a=(alpha/a[2])/(alpha/a[2]+(1-alpha)/b[2]);raw_b=1-raw_a
pacc=raw_a*a[2]+raw_b*b[2]
facc=raw_a*(1-a[0])*a[2]+raw_b*(1-b[0])*b[2]
assert facc==(1-eps)*pacc
m_a=1 if a[3]=='1' else 2 if a[3].startswith('2') else 3
m_b=1 if b[3]=='1' else 2 if b[3].startswith('2') else 3
assert (raw_a*m_a+raw_b*m_b)/pacc==co
assert raw_a+raw_b==1 and raw_a>=0 and raw_b>=0
# Exact residual CP decomposition: every native corrected map is a mixture
# of I,U,Y,V, and one-call H correction gives (id+Ad_Z)/2. All
# nonidentity unitary Choi vectors are orthogonal to Phi. Therefore each
# normalized accepted mixture has half-diamond error exactly 1-facc/pacc.
baseline=1+2*(6+x)/(2*x-1)*(F(1,2)-eps)
saving=1-co/baseline
print('t,x,q,baseline',t,x,q,baseline,float(baseline))
print('best',co,float(co),'alpha',alpha,'vertices',a,b,'raw',raw_a,raw_b,'pacc',pacc,'saving',saving,float(saving))
(_output_dir/'physical_proportional_certificate.json').write_text(json.dumps({'t':str(t),'eta2':str(x),'q':str(q),'timings':'T1=tau,T2=2tau,T3=3tau','baseline12':str(baseline),'epsilon':str(eps),'achieved_error':str(1-facc/pacc),'Bell_numerator':str(facc),'timecost_in_tau_units':str(co),'expected_queries_per_raw_trial':str(raw_a*m_a+raw_b*m_b),'codewords':[format(k,'03b')for k in CODES],'flag_p_f':[[str(z)for z in w]for w in weights],'optimum_within_native_subset_hull':str(co),'policy_vertices':[a[3],b[3]],'accepted_weight_first':str(alpha),'raw_weight_first':str(raw_a),'raw_weight_second':str(raw_b),'p_acc':str(pacc),'relative_saving_vs_optimal12':str(saving),'status':'EXACT CONSTRUCTIVE COUNTEREXAMPLE; NOT A GLOBAL OPTIMUM CLAIM'},indent=2))
