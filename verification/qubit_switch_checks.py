"""Finite-matrix verification of actual eight-order qubit switch branches.

Each noise history gets one full block-controlled switch Kraus W_a. Preparation
and heralding are applied as matrix contractions, and the complete prepared-switch
Choi output (including the unmeasured controller) is retained for the gauge test.
Uses only NumPy; no centrality substitution is made in constructing the channel.
"""
from itertools import product, permutations
from pathlib import Path
import argparse
import json
import numpy as np
from depolarizing_probability import QUBIT_TERMS

I=np.eye(2,dtype=complex)
PAULI=np.array([I,[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],dtype=complex)
WORDS=[tuple(x-1 for x in word) for word,_ in QUBIT_TERMS]
COEFFICIENTS=np.array([c for _,c in QUBIT_TERMS],dtype=complex)


def controllers(coefficients):
    l1=np.sum(np.abs(coefficients))
    alpha=np.sqrt(np.abs(coefficients)/l1).astype(complex)
    beta=np.exp(-1j*np.angle(coefficients))*alpha
    return alpha,beta,l1


def choi(kraus):
    vectors=np.array(kraus).reshape(len(kraus),-1)
    return vectors.T@vectors.conj()


def random_channel(rng,rank=4):
    z=rng.normal(size=(2*rank,2))+1j*rng.normal(size=(2*rank,2))
    q,_=np.linalg.qr(z)
    return q.reshape(rank,2,2)


def gauge_change(kraus,rng,new_size):
    z=rng.normal(size=(new_size,len(kraus)))+1j*rng.normal(size=(new_size,len(kraus)))
    v,_=np.linalg.qr(z)
    return np.einsum('ba,aij->bij',v,kraus)


def switch_kraus(slots,words=WORDS):
    control_size=len(words)
    for indices in product(*(range(len(slot)) for slot in slots)):
        w=np.zeros((2*control_size,2*control_size),complex)
        selected=[slot[a] for slot,a in zip(slots,indices)]
        for k,word in enumerate(words):
            matrix=I.copy()
            for j in word: matrix=matrix@selected[j]
            w[np.ix_([k,control_size+k],[k,control_size+k])]=matrix
        yield w,selected


def covariance(kraus):
    # K_a = sum_mu k_(a,mu) sigma_mu, with tr(sigma_mu sigma_nu)=2 delta_mu,nu.
    vectors=np.einsum('ixy,ayx->ai',PAULI[1:],np.array(kraus))/2
    return vectors.T@vectors.conj()


EPS=np.zeros((3,3,3),int)
for perm in permutations(range(3)):
    EPS[perm]=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))


def mixed_cross_covariance(q1,q2):
    return np.einsum('ikl,jmn,km,ln->ij',EPS,EPS,q1,q2)


def analytic_probability(slots):
    q=[covariance(k) for k in slots]
    s12=mixed_cross_covariance(q[0],q[1])
    s34=mixed_cross_covariance(q[2],q[3])
    return float(np.einsum('ij,ij->',s12,s34).real)


def adjugate3(matrix):
    cofactor=np.empty((3,3),complex)
    for i in range(3):
        for j in range(3):
            minor=np.delete(np.delete(matrix,i,axis=0),j,axis=1)
            cofactor[i,j]=(-1)**(i+j)*np.linalg.det(minor)
    return cofactor.T


def evaluate(slots,reference_state):
    alpha,beta,l1=controllers(COEFFICIENTS)
    s=len(alpha)
    prepare=np.kron(I,alpha[:,None])
    herald=np.kron(I,beta.conj()[None,:])
    omega=I.reshape(-1)
    ideal_choi=np.outer(omega,omega.conj())
    # Unnormalized maximally entangled target-reference input, with independent control.
    input_vector=np.einsum('sr,c->scr',I,alpha).reshape(2*s,2)
    joint=np.zeros((4*s,4*s),complex)
    success_choi=np.zeros((4,4),complex)
    tp=np.zeros((2*s,2*s),complex)
    success_effect=np.zeros((2,2),complex)
    output_reference=np.zeros_like(reference_state)
    reference_dim=reference_state.shape[0]//2
    scalar_residual=0.; count=0
    for w,_ in switch_kraus(slots):
        count+=1
        tp+=w.conj().T@w
        b=herald@w@prepare
        v=b.reshape(-1)
        success_choi+=np.outer(v,v.conj())
        scalar_residual=max(scalar_residual,float(np.linalg.norm(b-np.trace(b)*I/2)))
        success_effect+=b.conj().T@b
        full_output=(w@input_vector).reshape(-1)
        joint+=np.outer(full_output,full_output.conj())
        br=np.kron(b,np.eye(reference_dim))
        output_reference+=br@reference_state@br.conj().T
    p=float(np.trace(success_choi).real/2)
    postselect_reference=np.kron(np.kron(I,beta.conj()[None,:]),I)
    choi_from_joint=postselect_reference@joint@postselect_reference.conj().T
    eta=[max(0.,2*float(np.linalg.eigvalsh(choi(k))[0])) for k in slots]
    bound=(3/64)*float(np.prod(eta))
    stats={
        'noise_histories':count,'success_probability':p,
        'analytic_probability':analytic_probability(slots),
        'cp_noise_floor_parameters':eta,'cp_noise_floor_bound':bound,
        'switch_tp_error':float(np.linalg.norm(tp-np.eye(2*s))),
        'maximum_scalar_kraus_residual':scalar_residual,
        'identity_choi_error':float(np.linalg.norm(success_choi-p*ideal_choi)),
        'choi_joint_projection_error':float(np.linalg.norm(success_choi-choi_from_joint)),
        'input_independence_effect_error':float(np.linalg.norm(success_effect-p*I)),
        'reference_state_error':float(np.linalg.norm(output_reference-p*reference_state)),
        'complete_prepared_switch_choi_trace':float(np.trace(joint).real),
        'controller_dimension':s,'controller_norms':[float(np.linalg.norm(alpha)),float(np.linalg.norm(beta))],
    }
    for key in ['switch_tp_error','maximum_scalar_kraus_residual','identity_choi_error',
                'choi_joint_projection_error','input_independence_effect_error','reference_state_error']:
        assert stats[key]<1e-11,(key,stats[key])
    assert abs(p-stats['analytic_probability'])<1e-12
    assert p+1e-13>=bound
    assert abs(np.trace(joint).real-2)<1e-11
    return stats,joint,success_choi


def verify():
    rng=np.random.default_rng(20261028)
    z=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
    rho=z@z.conj().T; rho/=np.trace(rho)
    dep=PAULI/2
    pauli=np.sqrt(np.array([.17,.31,.23,.29]))[:,None,None]*PAULI
    gamma=.43
    ad=np.array([[[1,0],[0,np.sqrt(1-gamma)]],[[0,np.sqrt(gamma)],[0,0]]],complex)
    generic=random_channel(rng,4)
    distinct=[random_channel(rng,4) for _ in range(4)]
    cases={
        'completely_depolarizing':[dep]*4,
        'pauli_full_rank':[pauli]*4,
        'generic_full_rank':[generic]*4,
        'four_distinct_full_rank':distinct,
        'amplitude_damping_boundary':[ad]*4,
        'unitary_boundary':[np.array([PAULI[2]])]*4,
    }
    results={}; cached={}
    for name,slots in cases.items():
        stats,joint,success=evaluate(slots,rho)
        if all(k is slots[0] for k in slots):
            adj=adjugate3(covariance(slots[0]))
            p_adj=float((4*np.trace(adj@adj.T)).real)
            stats['adjugate_formula_probability']=p_adj
            assert abs(p_adj-stats['success_probability'])<1e-12
        results[name]=stats; cached[name]=(joint,success)
    assert abs(results['completely_depolarizing']['success_probability']-3/64)<1e-14
    assert results['amplitude_damping_boundary']['success_probability']<1e-26
    assert results['unitary_boundary']['success_probability']<1e-26

    # Independent rectangular Kraus isometries in every slot, including extra Kraus operators.
    gauged=[gauge_change(k,rng,6) for k in distinct]
    stats,joint,success=evaluate(gauged,rho)
    old_joint,old_success=cached['four_distinct_full_rank']
    stats['slot_choi_gauge_errors']=[float(np.linalg.norm(choi(k)-choi(g))) for k,g in zip(distinct,gauged)]
    stats['complete_prepared_switch_gauge_error']=float(np.linalg.norm(joint-old_joint))
    stats['success_choi_gauge_error']=float(np.linalg.norm(success-old_success))
    assert max(stats['slot_choi_gauge_errors'])<1e-12
    assert stats['complete_prepared_switch_gauge_error']<1e-12
    assert stats['success_choi_gauge_error']<1e-12
    results['rectangular_kraus_gauge_change']=stats

    # Phase convention check for a noncentral polynomial with unequal complex coefficients.
    coeff=rng.normal(size=8)+1j*rng.normal(size=8)
    alpha,beta,l1=controllers(coeff)
    first_w,selected=next(switch_kraus(distinct))
    branch=np.kron(I,beta.conj()[None,:])@first_w@np.kron(I,alpha[:,None])
    polynomial=np.zeros((2,2),complex)
    for c,word in zip(coeff,WORDS):
        term=I.copy()
        for j in word: term=term@selected[j]
        polynomial+=c*term
    phase_error=float(np.linalg.norm(branch-polynomial/l1))
    assert phase_error<1e-14
    results['coefficient_phase_check']={'l1_norm':float(l1),'branch_error':phase_error,
                                       'sum_abs_alpha_beta':float(np.sum(np.abs(alpha*beta.conj())))}
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Optional JSON report path')
    args=parser.parse_args()
    results=verify()
    if args.output:
        args.output.write_text(json.dumps(results,indent=2)+'\n')
    summary={name:{k:v for k,v in row.items() if k in ('success_probability','cp_noise_floor_bound',
              'identity_choi_error','complete_prepared_switch_gauge_error','branch_error')}
             for name,row in results.items()}
    print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
