from fractions import Fraction as F
from pathlib import Path
import itertools as it,json,math
import numpy as np
import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent,
                     help='Directory for regenerated exact certificate records')
_args=_parser.parse_args() if __name__=='__main__' else _parser.parse_args([])
ROOT=_args.output_dir.resolve()
ROOT.mkdir(parents=True,exist_ok=True)
ORDERS=list(it.permutations(range(3)))
CODES=[0,1,4,6,3,7]
ALT_CODES=[0,3,5,2,4,1]
ID=np.eye(2,dtype=int)
KINT=[np.array([[1,0],[1,0]],dtype=int),np.array([[0,1],[0,-1]],dtype=int)]
IVECS=[]
for hist in it.product(range(2),repeat=3):
 row=[]
 for order in ORDERS:
  a=ID
  for j in order:a=a@KINT[hist[j]]
  row.append(a.reshape(4))
 IVECS.append(np.array(row,dtype=int))
PAR=np.array([(-1)**sum(o[i]>o[j] for i in range(3) for j in range(i+1,3)) for o in ORDERS])
R=np.array([[int(a==tuple(reversed(b))) for b in ORDERS] for a in ORDERS])
T=np.full((6,6),F(1,6));S=np.array([[F(int(x),6) for x in rr] for rr in np.outer(PAR,PAR)])
O=np.array([[F(int(x),2) for x in rr] for rr in np.eye(6,dtype=int)-R])-S
E=np.array([[F(int(x),2) for x in rr] for rr in np.eye(6,dtype=int)+R])-T
PROJS=[O,E,T,S]
# Real representatives up to irrelevant global phase: I,U,Y,V
QS=[(ID,F(1)),(np.array([[-1,1],[1,1]]),F(1,2)),(np.array([[0,-1],[1,0]]),F(1)),(np.array([[1,1],[1,-1]]),F(1,2))]
BVS=[([1,0,0,1],F(1,2)),([-1,1,1,1],F(1,4)),([0,-1,1,0],F(1,2)),([1,1,1,-1],F(1,4))]

def coeffs(codes=CODES):
 d=np.array([[(a^b).bit_count() for b in codes] for a in codes]);ret=[]
 for P in PROJS:
  js=[]
  for h in range(4):
   J=np.zeros((4,4),dtype=object)
   for vs in IVECS:
    for i,j in it.product(range(6),repeat=2):
     if d[i,j]==h:J+=np.outer(vs[i],vs[j])*P[i,j]/96
   js.append(J)
  row=[]
  for Q,scale in QS:
   Q4=np.kron(Q,ID);cj=[scale*(Q4@j@Q4.T) for j in js]
   weights=[[sum(F(v[i]*v[j])*J[i,j]*s for i,j in it.product(range(4),repeat=2)) for J in cj] for v,s in BVS]
   # Certify each cubic coefficient is Bell diagonal after correction.
   for power,J in enumerate(cj):
    recon=np.zeros((4,4),dtype=object)
    for (v,scl),poly in zip(BVS,weights):
     recon+=np.outer(v,v)*scl*poly[power]
    assert np.array_equal(J,recon)
    pp=sum(poly[power] for poly in weights)
    assert np.array_equal(J[:2,:2]+J[2:,2:],np.eye(2,dtype=object)*pp/2)
   row.append(weights)
  ret.append(row)
 return ret

def evalpoly(poly,q):return sum(c*q**i for i,c in enumerate(poly))

def flagtable(q,codes=CODES):
 co=coeffs(codes);ps=[];fs=[];cor=[];diag=[]
 for data in co:
  r=max(range(4),key=lambda j:evalpoly(data[j][0],q))
  w=[evalpoly(poly,q) for poly in data[r]]
  ps.append(sum(w));fs.append(w[0]);cor.append(r);diag.append(w)
 return ps,fs,cor,diag

if __name__=='__main__':
 ret=coeffs()
 labels=['O','E','T','S'];out={}
 for label,row in zip(labels,ret):
  print(label)
  out[label]={name:[[str(c) for c in poly] for poly in weights] for name,weights in zip(['I','U','Y','V'],row)}
  for name,weights in zip(['I','U','Y','V'],row):
   print(name,'fpoly',list(map(str,weights[0])),'ppoly',list(map(str,np.sum(weights,axis=0))))
 (ROOT/'pairwise_exact_flags.json').write_text(json.dumps(out,indent=2))
 ps,fs,cor,diag=flagtable(F(9,10));print('q.9',ps,fs,cor,diag,'score',sum(max(F(0),p+F(69,4)*(f-p/2)) for p,f in zip(ps,fs)))
