#!/usr/bin/env python3
"""Exact finite closure for m=1,2,3 quantum-order recovery.

Dependencies: NumPy, SciPy and Matplotlib, pinned in requirements.txt.
Use --output-dir DIR to keep regenerated JSON and figures separate from the
committed records; --skip-figures runs the complete checks without plotting.
The model and analytical converse are in the repository science pages.
"""
from dataclasses import dataclass,replace
from fractions import Fraction as F
from pathlib import Path
import json
import argparse
import itertools as it
import numpy as np
from scipy.optimize import linprog,brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import prior_three_call_certificate as prior

HERE=Path(__file__).resolve().parent

def flags2(eta):
    return [('minus',(2-eta)/4,(1+eta)/8),('plus',(2+eta)/4,(1+eta)/8)]

def flags3(eta):
    return [('O',(4-eta)/12,(1+2*eta)/12),('E',(4-eta)/12,(1+eta)/12),
            ('T',(1+2*eta)/6,(1+5*eta)/24),('S',(1-eta)/6,(1-eta)/24)]

def three_fidelity_numerator(p,eta):
    remain=p;out=0
    for _,prob,f in sorted([z for z in flags3(eta) if z[1]>0],key=lambda z:z[2]/z[1],reverse=True):
        take=min(remain,prob);out+=take*f/prob;remain-=take
        if remain<=0:break
    return out

@dataclass
class Policy:
    error:object
    cost:object
    probability:object
    calls:int
    name:str
    flags:tuple

def architecture_vertices(m,eta):
    if m==1:return [Policy(F(1,2),F(1),F(1),1,'1:H',('one',))]
    flags=flags2(eta) if m==2 else flags3(eta)
    flags=sorted([z for z in flags if z[1]>0],key=lambda z:z[2]/z[1],reverse=True)
    prob=f=0;out=[];names=[]
    for name,p,g in flags:
        prob+=p;f+=g;names.append(name)
        out.append(Policy(1-f/prob,m/prob,prob,m,f'{m}:'+''.join(names),tuple(names)))
    return out

def vertices(eta2,eta3,include3=True,cost_weights=None):
    out=architecture_vertices(1,0)+architecture_vertices(2,eta2)+(architecture_vertices(3,eta3) if include3 else [])
    if cost_weights is not None:out=[replace(p,cost=cost_weights[p.calls-1]/p.probability) for p in out]
    return out

def lower_hull(vs):
    # Pareto prune, coalesce equal errors, then take increasing-slope lower hull.
    monotone=[]
    for p in sorted(vs,key=lambda p:(p.error,p.cost)):
        if monotone and p.error==monotone[-1].error:continue
        if monotone and p.cost>=monotone[-1].cost:continue
        monotone.append(p)
    hull=[]
    for p in monotone:
        while len(hull)>=2:
            a,b=hull[-2:]
            cross=(b.error-a.error)*(p.cost-b.cost)-(b.cost-a.cost)*(p.error-b.error)
            if cross<=0:hull.pop()
            else:break
        hull.append(p)
    return hull

def optimal_cost(epsilon,eta2,eta3,include3=True,return_plan=False,cost_weights=None):
    hs=lower_hull(vertices(eta2,eta3,include3,cost_weights))
    if epsilon<hs[0].error:return (None,[]) if return_plan else None
    if epsilon>=hs[-1].error:
        return (hs[-1].cost,[(hs[-1],1)]) if return_plan else hs[-1].cost
    for a,b in zip(hs[:-1],hs[1:]):
        if a.error<=epsilon<=b.error:
            alpha=(b.error-epsilon)/(b.error-a.error)
            cost=alpha*a.cost+(1-alpha)*b.cost
            return (cost,[(a,alpha),(b,1-alpha)]) if return_plan else cost
    raise RuntimeError('hull interpolation failed')

def trial_plan(plan):
    norm=sum(alpha/p.probability for p,alpha in plan)
    return [(p,(alpha/p.probability)/norm) for p,alpha in plan],1/norm

def flow_lp(epsilon,eta2,eta3,cost_weights=(1,2,3)):
    # Charnes-Cooper exact finite LP: variables normalized accepted flag flows
    # u_(m,z), and v_m=raw architecture probability / aggregate acceptance.
    all_flags=[(1,'one',1.,.5)]
    all_flags += [(2,n,float(p),float(f)) for n,p,f in flags2(eta2)]
    all_flags += [(3,n,float(p),float(f)) for n,p,f in flags3(eta3)]
    n=len(all_flags);objective=np.zeros(n+3);objective[n:]=cost_weights
    caps=[]
    for i,(m,_,_,_) in enumerate(all_flags):
        row=np.zeros(n+3);row[i]=1;row[n+m-1]=-1;caps.append(row)
    quality=np.zeros(n+3);quality[:n]=[-float(f) for _,_,_,f in all_flags]
    ps=np.zeros(n+3);ps[:n]=[float(p) for _,_,p,_ in all_flags]
    return linprog(objective,A_ub=np.vstack(caps+[quality]),b_ub=np.r_[np.zeros(n),-(1-epsilon)],
                   A_eq=[ps],b_eq=[1.],bounds=[(0,None)]*(n+3),method='highs')

def baseline_cost(epsilon,eta2):
    if epsilon>=.5:return 1.
    if eta2<=.5:return None
    e0=3*(1-eta2)/(2*(2-eta2))
    if epsilon<e0:return None
    return 1+2*(6+eta2)*(.5-epsilon)/(2*eta2-1)

def cost_phase_score(x,z):
    if x<=.5:return None
    b=2*(6+x)/(2*x-1)
    ks=[(8-2*z+b*(5*z-2))/24,(8-2*z+b*(3*z-2))/24,(4+8*z+b*(z-1))/24]
    return sum(max(v,0) for v in ks)-3

def any_three_activation(x,z):
    if x<=.5:return z>.4
    return z>min(2*x/(3-x),(57*x+4)/(13*x+52))

def priority_cost(z):
    if z<=F(79,91):return F(101,32)
    if z<F(15,17):return (789-680*z)/(135-83*z)
    return (17*z+3)/(18*z-10)

def policy_choi(policy,eta2,eta3):
    if policy.calls==1:
        vz=np.kron(prior.Z,prior.I)@prior.phi
        return .5*prior.PHI+.5*np.outer(vz,vz.conj())
    if policy.calls==2:
        pm,pp=(2-eta2)/4,(2+eta2)/4;ell,s=(1+eta2)/8,(1-eta2)/8
        jm=prior.BELL@np.diag([ell,s,s,s])@prior.BELL.conj().T
        jp=prior.BELL@np.diag([ell,s,ell,ell])@prior.BELL.conj().T
        j=jm if len(policy.flags)==1 else jm+jp
        return j/float(policy.probability)
    j=np.zeros((4,4),complex)
    for name,p,c in zip(['O','E','T','S'],prior.PROJECTORS,prior.CORRECTIONS):
        if name in policy.flags:j+=prior.branch_choi(p,float(eta3),c)
    return j/float(policy.probability)

def generate_figures(output_dir=HERE):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(11,3.9),layout='constrained')
    es=np.linspace(.075,.5,500)
    for include,label,color in [(False,'Optimal one/two-call mixture','#b34f37'),(True,'Optimal one/two/three-call mixture','#206f8a')]:
        ys=[optimal_cost(e,.9,.9,include) for e in es]
        ax[0].plot(es,[np.nan if c is None else c for c in ys],label=label,color=color,lw=2)
    ax[0].scatter([.375],[183/62],color='#206f8a',s=28,zorder=4)
    ax[0].set(xlabel='Half-diamond error budget ε',ylabel='Expected noisy calls / accepted trial',title='Exact convex closure at η₂ = η₃ = 0.9',ylim=(.7,13))
    ax[0].legend(frameon=False,fontsize=8.5)
    zs=np.linspace(.82,1,500)
    ax[1].plot(zs,[optimal_cost(.375,.9,z) for z in zs],color='#206f8a',lw=2,label='Optimal ≤3-call mixture')
    ax[1].axhline(101/32,color='#b34f37',ls='--',label='Optimal ≤2-call baseline')
    ax[1].axvline(79/91,color='#89929a',ls=':',lw=1)
    ax[1].set(xlabel='Three-call coherence η₃',ylabel='Expected noisy calls / accepted trial',title='Unequal noise: η₂ = 0.9, ε ≤ 3/8',ylim=(2.4,3.25))
    ax[1].legend(frameon=False,fontsize=8.5)
    fig.savefig(output_dir/'closure_cost_and_noise.png',dpi=210,facecolor='white')
    fig2,ax=plt.subplots(1,2,figsize=(11,3.9),layout='constrained')
    xs=np.linspace(.501,1,400);zs=np.linspace(.35,1,380)
    xx,zz=np.meshgrid(xs,zs);score=np.array([[cost_phase_score(x,z) for x in xs] for z in zs])
    ax[0].contourf(xx,zz,score>0,levels=[-.5,.5,1.5],colors=['#f1ede8','#b9dce7'])
    ax[0].contour(xx,zz,score,levels=[0],colors=['#206f8a'],linewidths=1.8)
    ax[0].plot(xs,2*xs/(3-xs),color='#65518a',ls='--',label='Accuracy-extension boundary')
    ax[0].plot(xs,xs,color='#6d6d6d',ls=':',label='Equal coherence')
    ax[0].plot(xs,xs**1.5,color='#b34f37',label='Call-proportional duration')
    ax[0].set(xlabel='Two-call coherence η₂',ylabel='Three-call coherence η₃',title='Exact common-error call-cost advantage',ylim=(.35,1))
    ax[0].legend(frameon=False,fontsize=8,loc='lower right')
    rat=np.linspace(1,1.65,350)
    ax[1].plot(rat,[optimal_cost(.375,.9,.9**a) for a in rat],color='#206f8a',lw=2)
    ax[1].axhline(101/32,color='#b34f37',ls='--')
    amax=np.log(79/91)/np.log(.9)
    ax[1].axvline(amax,color='#89929a',ls=':',lw=1)
    ax[1].text(1.03,3.12,'Conditional scalar dephasing: ηₘ = exp(−γTₘ)',fontsize=8.5)
    ax[1].set(xlabel='Duration ratio T₃ / T₂',ylabel='Optimal calls / accepted trial',title='Duration budget at η₂ = 0.9, ε ≤ 3/8',ylim=(2.8,3.22))
    fig2.savefig(output_dir/'closure_advantage_region.png',dpi=210,facecolor='white')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=HERE,help='directory for regenerated results and figures (default: beside this script)')
    parser.add_argument('--skip-figures',action='store_true',help='run every check without generating plots')
    args=parser.parse_args()
    output_dir=args.output_dir.resolve();output_dir.mkdir(parents=True,exist_ok=True)
    prior.exact_certificate()
    max_choi=max_lp=max_hull=0.;count=0
    for eta in [0,.1,1/7,.4,.5,7/13,.6,.9,1]:
        for p in np.linspace(.005,1,61):
            fs=flags3(eta)
            lp=linprog([-v[2] for v in fs],A_eq=[[v[1] for v in fs]],b_eq=[p],bounds=[(0,1)]*4,method='highs')
            assert lp.success
            f=three_fidelity_numerator(p,eta);max_lp=max(max_lp,abs(f+lp.fun))
            j=sum(a*prior.branch_choi(P,eta,c) for a,P,c in zip(lp.x,prior.PROJECTORS,prior.CORRECTIONS))
            np.testing.assert_allclose(np.trace(j),p,atol=1e-12)
            e=np.abs(np.linalg.eigvalsh(j/p-prior.PHI)).sum()/2
            max_choi=max(max_choi,abs(e-(1-f/p)))
            np.testing.assert_allclose(np.einsum('iaib->ab',j.reshape(2,2,2,2)),p*prior.I/2,atol=2e-13)
            count+=1
    rng=np.random.default_rng(29102026)
    for _ in range(350):
        x,z,eps=rng.uniform(),rng.uniform(),rng.uniform(0,.7)
        expected=optimal_cost(eps,x,z);lp=flow_lp(eps,x,z)
        if expected is None:assert lp.status==2
        else:
            assert lp.success
            max_hull=max(max_hull,abs(lp.fun-expected))
            c,plan=optimal_cost(eps,x,z,return_plan=True)
            j=sum(float(a)*policy_choi(pol,x,z) for pol,a in plan)
            actual=np.abs(np.linalg.eigvalsh(j-prior.PHI)).sum()/2
            assert actual<=eps+2e-11
    for _ in range(100):
        x,z,eps=rng.uniform(),rng.uniform(),rng.uniform(0,.7)
        weights=rng.uniform(.1,5,size=3)
        expected=optimal_cost(eps,x,z,cost_weights=weights);lp=flow_lp(eps,x,z,weights)
        if expected is None:assert lp.status==2
        else:
            assert lp.success
            max_hull=max(max_hull,abs(lp.fun-expected))
    # Arbitrary complex accepted controller Grams and continuum corrections.
    max_global_violation=0.
    for _ in range(250):
        w=rng.dirichlet(np.ones(6));effects=[]
        for y in range(4):
            a=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6));effects.append(a@a.conj().T)
        d,v=np.linalg.eigh(sum(effects));inv=(v*(1/np.sqrt(d)))@v.conj().T
        effects=[inv@a@inv for a in effects];eta=rng.uniform();j=np.zeros((4,4),complex)
        for f in effects[:3]:
            g=np.sqrt(w)[:,None]*f.T*np.sqrt(w)[None,:]
            raw=prior.choi(eta*g+(1-eta)*np.diag(np.diag(g)))
            h=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2));q,_=np.linalg.qr(h);q=np.kron(q,prior.I)
            j+=q@raw@q.conj().T
        p=np.trace(j).real;f=np.vdot(prior.phi,j@prior.phi).real
        max_global_violation=max(max_global_violation,f-three_fidelity_numerator(p,eta))
    assert max_global_violation<1e-12
    # Exact rational priority hull, trial policy, and complete cost formula.
    x=F(9,10);eps=F(3,8)
    hs=lower_hull(vertices(x,x));assert [p.error for p in hs]==[F(3,31),F(3,22),F(87,236),F(1,2)]
    c,plan=optimal_cost(eps,x,x,return_plan=True);assert c==F(183,62)
    trials,pacc=trial_plan(plan);assert pacc==F(62,63)
    trialdict={p.calls:t for p,t in trials};assert trialdict=={3:F(20,21),1:F(1,21)}
    baseline=optimal_cost(eps,x,x,False);assert baseline==F(101,32)
    saving=1-c/baseline;assert saving==F(203,3131)
    for z in [F(i,200) for i in range(201)]+[F(79,91),F(15,17),F(9,10)]:
        assert optimal_cost(eps,x,z)==priority_cost(z)
    # Exact phase formula vs cost LP on a useful common-error slice.
    phase_count=0
    for x in np.linspace(.51,1,45):
        e0=3*(1-x)/(2*(2-x));e=(e0+.5)/2
        for z in np.linspace(0,1,55):
            score=cost_phase_score(x,z);c3=optimal_cost(e,x,z);c12=baseline_cost(e,x)
            if abs(score)<1e-8:continue
            got=c3 is not None and c3<c12-1e-9
            assert got==(score>0),(x,z,score,c3,c12)
            phase_count+=1
    for _ in range(1500):
        x,z=rng.uniform(),rng.uniform();vs3=architecture_vertices(3,z)
        floor=min(float(p.error) for p in vertices(x,z,False))
        quality=min(float(p.error) for p in vs3)<floor-1e-12
        cost=(x>.5 and cost_phase_score(x,z)>1e-10)
        assert any_three_activation(x,z)==(quality or cost)
    eta_low=(35-np.sqrt(985))/6;eta_high=(5+np.sqrt(233))/26
    for eta in np.linspace(.5001,1,1000):
        assert (cost_phase_score(eta,eta)>0)==(eta<eta_low or eta>eta_high)
    troot=brentq(lambda t:13*t**5+52*t**3-57*t*t-4,.95,1,xtol=5e-15)
    xr=troot*troot;aratio=np.log(79/91)/np.log(.9)
    physical_poly=lambda q:169*q**5+1352*q**4+2704*q**3-3249*q*q-456*q-16
    physical_lo=F('0.9199910805350899');physical_hi=F('0.9199910805350900')
    assert physical_poly(physical_lo)<0<physical_poly(physical_hi)
    assert abs(cost_phase_score(xr,xr**1.5))<1e-10
    np.testing.assert_allclose(optimal_cost(.375,.9,.9**1.5),baseline_cost(.375,.9),atol=1e-13)
    assert max_lp<1e-10 and max_choi<1e-10 and max_hull<1e-8
    output={'status':'PASS','three_call_direct_grid_points':count,'raw_flow_lp_hull_cases':350,'additional_resource_weight_cases':100,
            'arbitrary_complex_gram_strategies':250,'phase_grid_cases':phase_count,
            'activation_random_cases':1500,'max_three_knapsack_lp_residual':max_lp,
            'max_direct_choi_error_residual':max_choi,'max_raw_flow_vs_hull_residual':max_hull,
            'positive_global_witness_violation':max(0.,max_global_violation),
            'priority':{'eta2':'9/10','eta3':'9/10','epsilon':'3/8','optimal_cost':str(c),
                        'baseline_cost':str(baseline),'relative_saving':str(saving),'acceptance':str(pacc),
                        'raw_trial_probabilities':{str(k):str(v) for k,v in trialdict.items()}},
            'equal_noise_reentrant_thresholds':[eta_low,eta_high],
            'unequal_noise_priority_threshold':'79/91','priority_duration_ratio_limit':aratio,
            'per_call_duration_ratio_1_5_coherence_threshold':xr,
            'per_call_exponent_threshold':-np.log(troot),
            'exact_physical_coherence_root_bracket':[str(physical_lo),str(physical_hi)],
            'conditional_physics_warning':'Scalar logical-order dephasing model only; local physical-qubit errors need a different channel.'}
    (output_dir/'complete_closure_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
    if not args.skip_figures:generate_figures(output_dir)

if __name__=='__main__':main()
