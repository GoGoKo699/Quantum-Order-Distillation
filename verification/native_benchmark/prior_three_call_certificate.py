#!/usr/bin/env python3
"""Exact/numerical checks of the three-call advantage. Requires numpy, scipy, matplotlib.
Imported by verify_complete_closure.py for native projectors, branch Choi maps
and exact rational certificates. Direct execution supports --output-dir DIR.
Words are matrix-product order; Choi states use a normalized Bell input.
"""
from pathlib import Path
from fractions import Fraction as F
import itertools as it
import json
import argparse
import numpy as np
from scipy.optimize import linprog
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.])
U=(X-Z)/np.sqrt(2);V=(X+Z)/np.sqrt(2);Q=[I,U,Y,V]
K=[np.array([[1,0],[1,0]],complex)/np.sqrt(2),np.array([[0,1],[0,-1]],complex)/np.sqrt(2)]
ORDERS=list(it.permutations(range(3))); r=6
parity=np.array([(-1)**sum(o[i]>o[j] for i in range(3) for j in range(i+1,3)) for o in ORDERS])
REV=np.array([[int(a==tuple(reversed(b))) for b in ORDERS] for a in ORDERS])
T=np.ones((6,6))/6;S=np.outer(parity,parity)/6
O=(np.eye(6)-REV)/2-S;E=(np.eye(6)+REV)/2-T
PROJECTORS=[O,E,T,S];RANKS=[2,2,1,1];CORRECTIONS=[U,Y,I,I]
phi=np.array([1,0,0,1],complex)/np.sqrt(2);PHI=np.outer(phi,phi.conj());JDEP=np.eye(4)/4
BELL=np.stack([np.kron(q,I)@phi for q in Q],axis=1)

def words():
    out=[]
    for h in it.product(range(2),repeat=3):
        row=[]
        for o in ORDERS:
            a=I.copy()
            for j in o:a=a@K[h[j]]
            row.append(a)
        out.append(row)
    return np.array(out)
A=words()

def choi(g):
    vectors=A.reshape(8,6,4)/np.sqrt(2)
    return np.einsum('ij,hia,hjb->ab',g,vectors,vectors.conj())

def branch_choi(projector,eta,correction=I):
    g=(eta*projector+(1-eta)*np.diag(np.diag(projector)))/6
    j=choi(g);q=np.kron(correction,I)
    return q@j@q.conj().T

def two_pmax(epsilon,eta):
    e0=3*(1-eta)/(2*(2-eta));ed=(3-eta)/4
    if epsilon<e0-1e-13:return 0.
    if epsilon>=ed-1e-13:return 1.
    return eta*(1+eta)/(2*(3+eta-2*(2+eta)*epsilon))

def three_floor(eta):
    return 3*(1+eta)/(4*(1+2*eta)) if eta<=1/7 else 3*(1-eta)/(4-eta)

def mixed_one_two_cost(epsilon,eta):
    if epsilon>=.5:return 1.
    e0=3*(1-eta)/(2*(2-eta))
    if eta<=.5 or epsilon<e0:return None
    return 1+2*(6+eta)*(.5-epsilon)/(2*eta-1)

def family(eta,epsilon):
    js=[branch_choi(p,eta,c) for p,c in zip(PROJECTORS,CORRECTIONS)]
    probs=np.array([np.trace(j).real for j in js]);bells=np.array([np.vdot(phi,j@phi).real for j in js])
    lp=linprog(-probs,A_ub=[(1-epsilon)*probs-bells],b_ub=[0.],bounds=[(0,1)]*4,method='highs')
    assert lp.success
    acc=sum(a*j for a,j in zip(lp.x,js));p=np.trace(acc).real
    if p<1e-10:return 0.,None,lp.x
    error=np.abs(np.linalg.eigvalsh(acc/p-PHI)).sum()/2
    assert error<=epsilon+1e-10
    return float(p),float(error),lp.x

def mm(a,b):return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def outer(v):return [[v[i]*v[j] for j in range(4)] for i in range(4)]
def exact_certificate():
    t=[[F(1,6)]*6 for _ in range(6)]
    s=[[F(int(parity[i]*parity[j]),6) for j in range(6)] for i in range(6)]
    o=[[F(int(i==j)-int(REV[i,j]),2)-s[i][j] for j in range(6)] for i in range(6)]
    e=[[F(int(i==j)+int(REV[i,j]),2)-t[i][j] for j in range(6)] for i in range(6)]
    m=[[[1,0],[1,0]],[[0,1],[0,-1]]]
    jbasis=[[[F(a,2) for a in row] for row in outer([1,0,0,1])],
            [[F(a,4) for a in row] for row in outer([-1,1,1,1])],
            [[F(a,2) for a in row] for row in outer([0,-1,1,0])],
            [[F(a,4) for a in row] for row in outer([1,1,1,-1])]]
    coefficients=[[0,F(1,4),0,0],[0,0,F(1,6),F(1,12)],[F(1,4),0,F(1,12),F(1,6)],[0,0,0,0]]
    # Integer word data give exact probability and Bell-witness Gram matrices.
    ivecs=[]
    for hist in it.product(range(2),repeat=3):
        row=[]
        for order in ORDERS:
            a=[[1,0],[0,1]]
            for k in order:a=mm(a,m[hist[k]])
            row.append([x for rr in a for x in rr])
        ivecs.append(row)
    he=[[F(sum(sum(v[i][a]*v[j][a] for a in range(4)) for v in ivecs),16) for j in range(6)] for i in range(6)]
    assert he==[[3*t[i][j]+F(3,4)*(o[i][j]+e[i][j]) for j in range(6)] for i in range(6)]
    ints=[[1,0,0,1],[-1,1,1,1],[0,-1,1,0],[1,1,1,-1]]
    traces=[[[sum(q[a]*v[i][a] for a in range(4)) for i in range(6)] for v in ivecs] for q in ints]
    wants=[[[F(3,2)*t[i][j] for j in range(6)] for i in range(6)],
           [[F(3,4)*o[i][j] for j in range(6)] for i in range(6)],
           [[F(1,2)*(t[i][j]+e[i][j]) for j in range(6)] for i in range(6)],
           [[t[i][j]+F(1,4)*e[i][j] for j in range(6)] for i in range(6)]]
    for tr,denom,want in zip(traces,[32,64,32,64],wants):
        got=[[F(sum(v[i]*v[j] for v in tr),denom) for j in range(6)] for i in range(6)]
        assert got==want
    for a,b in it.combinations(range(4),2):
        if (a,b)==(1,3):continue
        assert all(sum(traces[a][h][i]*traces[b][h][j] for h in range(8))==0 for i,j in it.product(range(6),repeat=2))
    cert=[]
    for name,p,w in zip(['O','E','T','S'],[o,e,t,s],coefficients):
        j=[[F(0) for _ in range(4)] for _ in range(4)]
        for hist in it.product(range(2),repeat=3):
            vec=[]
            for order in ORDERS:
                a=[[1,0],[0,1]]
                for k in order:a=mm(a,m[hist[k]])
                vec.append([x for row in a for x in row])
            for a,b in it.product(range(4),repeat=2):
                j[a][b]+=sum(p[k][l]*vec[k][a]*vec[l][b] for k,l in it.product(range(6),repeat=2))/96
        expected=[[sum(w[k]*jbasis[k][a][b] for k in range(4)) for b in range(4)] for a in range(4)]
        assert j==expected
        cert.append({'projector':name,'ideal_choi':[[str(x) for x in row] for row in j]})
    return {'ideal_branches':cert,'success_gram':[[str(x) for x in row] for row in he],'bell_witness_identities':'M_I=3T/2; M_U=3O/4; M_Y=(T+E)/2; M_V=T+E/4','only_nonzero_off_diagonal_Bell_block':'U,V'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=HERE,help='directory for optional standalone results and plot')
    args=parser.parse_args()
    output_dir=args.output_dir.resolve();output_dir.mkdir(parents=True,exist_ok=True)
    certificate=exact_certificate()
    for p,rank in zip(PROJECTORS,RANKS):
        np.testing.assert_allclose(p@p,p,atol=1e-14);np.testing.assert_allclose(np.trace(p),rank,atol=1e-14)
    np.testing.assert_allclose(sum(PROJECTORS),np.eye(6),atol=1e-14)
    for i,j in it.combinations(range(4),2):np.testing.assert_allclose(PROJECTORS[i]@PROJECTORS[j],0,atol=1e-14)
    expected_coeffs=[[0,.25,0,0],[0,0,1/6,1/12],[.25,0,1/12,1/6],[0,0,0,0]]
    for p,w in zip(PROJECTORS,expected_coeffs):
        np.testing.assert_allclose(choi(p/6),BELL@np.diag(w)@BELL.conj().T,atol=2e-14)
    # Exact success Gram identity and Bell witness identities (also checked rationally above).
    h=np.einsum('hiab,hjab->ij',A.conj(),A)/2
    np.testing.assert_allclose(h,3*T+.75*(O+E),atol=2e-14)
    witnesses=[]
    for q,want in zip(Q,[1.5*T,.75*O,.5*(T+E),T+.25*E]):
        trace=np.einsum('ab,hiab->hi',q.conj(),A)/2
        witness=trace.conj().T@trace;witnesses.append(witness)
        np.testing.assert_allclose(witness,want,atol=2e-14)
    max_residual=max_witness_violation=0.
    for eta in np.linspace(0,1,101):
        branches=[branch_choi(p,eta,c) for p,c in zip(PROJECTORS,CORRECTIONS)]
        j=sum(branches)
        expected=eta*BELL@np.diag([2/3,1/12,1/12,1/6])@BELL.conj().T+(1-eta)*JDEP
        max_residual=max(max_residual,float(np.max(np.abs(j-expected))))
        np.testing.assert_allclose(j,expected,atol=2e-14)
        np.testing.assert_allclose(np.einsum('iaib->ab',j.reshape(2,2,2,2)),I/2,atol=2e-14)
        error=np.abs(np.linalg.eigvalsh(j-PHI)).sum()/2
        np.testing.assert_allclose(error,.75-5*eta/12,atol=2e-14)
        floor_projector= T if eta<=1/7 else O
        floor_correction=I if eta<=1/7 else U
        jf=branch_choi(floor_projector,eta,floor_correction);pf=np.trace(jf).real
        np.testing.assert_allclose(np.abs(np.linalg.eigvalsh(jf/pf-PHI)).sum()/2,three_floor(eta),atol=2e-14)
        # Every raw corrected Choi block uses a valid input-independent acceptance map here.
        for b in branches:
            prob=np.trace(b).real
            np.testing.assert_allclose(np.einsum('iaib->ab',b.reshape(2,2,2,2)),prob*I/2,atol=2e-14)
        if eta>=.2:
            eps=.75-5*eta/12;p2=two_pmax(eps,eta)
            np.testing.assert_allclose(p2,3*(1+eta)/(7+5*eta),atol=2e-13)
            assert 2/p2>3
    rng=np.random.default_rng(30102026)
    for _ in range(250):
        v=rng.normal(size=4);v/=np.linalg.norm(v)
        q=v[0]*I+1j*(v[1]*X+v[2]*Y+v[3]*Z)
        trace=np.einsum('ab,hiab->hi',q.conj(),A)/2;m=trace.conj().T@trace
        std=O+E
        np.testing.assert_allclose(T@m@std,0,atol=2e-14)
        np.testing.assert_allclose(m@S,0,atol=2e-14)
        max_witness_violation=max(max_witness_violation,float(np.linalg.eigvalsh(std@m@std-.75*std).max()),float(np.trace(T@m).real-1.5))
    assert max_residual<1e-13 and max_witness_violation<1e-13
    table=[]
    for eps in [.1,.2,.25,.3,.375]:
        p3,e3,acc=family(.9,eps);p2=two_pmax(eps,.9)
        table.append({'eta':.9,'error_budget':eps,'two_call_optimal_probability':p2,
                      'two_call_optimal_cost':2/p2 if p2 else None,
                      'three_call_feasible_probability':p3,'three_call_feasible_cost':3/p3 if p3 else None,
                      'three_call_actual_error':e3,'acceptance_fractions_O_E_T_S':acc.tolist()})
    result={'status':'PASS','exact_rational_branch_certificates':certificate,'noise_grid_points':101,
            'random_correction_witness_checks':250,'max_deterministic_choi_residual':max_residual,
            'max_positive_witness_bound_violation':max(0.,max_witness_violation),
            'cost_saving_eta_09':1-3/(230/57),
            'stronger_one_two_baseline_cost_eta_09':mixed_one_two_cost(.375,.9),
            'saving_over_stronger_baseline_eta_09':1-3/(101/32),
            'stronger_baseline_crossover_eta':(np.sqrt(129)-3)/10,'examples':table}
    # Accepted-ensemble and actual per-trial forms of the stronger baseline.
    jm=branch_choi(O,.9,U)
    # Use the two-call optimum's corrected map, not the three-call O herald.
    j2=(.9/4*PHI+.1/2*JDEP)/(.275)
    jz=.5*PHI+.5*np.outer(np.kron(Z,I)@phi,(np.kron(Z,I)@phi).conj())
    jmix=F(11,32)*j2+F(21,32)*jz
    jmix=np.asarray(jmix,dtype=complex)
    np.testing.assert_allclose(np.abs(np.linalg.eigvalsh(jmix-PHI)).sum()/2,3/8,atol=2e-14)
    assert F(21,61)+F(40,61)*F(11,40)==F(32,61)
    assert (F(21,61)+2*F(40,61))/F(32,61)==F(101,32)
    (output_dir/'three_call_advantage_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='exact_rational_branch_certificates'},indent=2))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(10.6,3.8),layout='constrained');etas=np.linspace(0,1,601)
    axs[0].plot(etas,3*(1-etas)/(2*(2-etas)),label='Two-call global minimum',color='#b34f37',lw=2)
    axs[0].plot(etas,[three_floor(e) for e in etas],label='Three-call global minimum',color='#206f8a',lw=2)
    axs[0].axvline(1/7,color='#89929a',ls=':',lw=1)
    axs[0].set(xlabel='Retained order coherence η',ylabel='Minimum conditional half-diamond error',title='A lower error becomes attainable',xlim=(0,1),ylim=(-.02,.78))
    axs[0].legend(frameon=False,fontsize=9)
    xs=np.linspace(.2,1,401);axs[1].plot(xs,2*(7+5*xs)/(3*(1+xs)),color='#b34f37',lw=2,label='Best two-call cost')
    mix=[mixed_one_two_cost(.75-5*e/12,e) for e in etas]
    axs[1].plot(etas,[np.nan if c is None else c for c in mix],color='#65518a',ls='--',lw=2,label='Best one/two-call mixture')
    axs[1].axhline(3,color='#206f8a',lw=2,label='Explicit three-call protocol')
    axs[1].axvspan(0,.2,alpha=.08,color='#206f8a');axs[1].text(.025,4.9,'Two calls\ninfeasible',fontsize=9)
    axs[1].set(xlabel='Retained order coherence η',ylabel='Noisy calls / accepted fresh trial',title='Matched error ε = 3/4 − 5η/12',xlim=(0,1),ylim=(.75,5.4))
    axs[1].legend(frameon=False,fontsize=9,loc='upper right')
    for ax in axs:ax.grid(alpha=.15)
    fig.savefig(output_dir/'three_call_advantage.png',dpi=210,facecolor='white')

if __name__=='__main__':main()
