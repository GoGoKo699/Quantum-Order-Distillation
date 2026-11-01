"""Run the exact and finite-matrix verification suite.

Default: Python standard library and NumPy only. Add --exhaustive to recompute
the 2880-word qutrit pair histogram with the supplied optional C++17 helper.
Existing coefficient exports and canonical result files are read, not rewritten.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

import check_access_model
import depolarizing_probability as depolarizing
import qubit_minimum
import qubit_switch_checks
import qutrit_degree8
import qutrit_trace_crosscheck
import regev_identity_checks
import regev_probabilities


HERE=Path(__file__).resolve().parent


def load_result(name):
    return json.loads((HERE/name).read_text())


def json_value(value):
    return json.loads(json.dumps(value))


def compare(actual,expected,location='result'):
    """Compare exact data exactly and floating records with explicit tolerance."""
    if isinstance(expected,float):
        if not isinstance(actual,(int,float)) or not math.isfinite(actual):
            raise AssertionError(f'{location}: expected a finite number')
        if not math.isclose(actual,expected,rel_tol=1e-10,abs_tol=1e-12):
            raise AssertionError(f'{location}: {actual!r} != {expected!r}')
    elif isinstance(expected,dict):
        if not isinstance(actual,dict) or actual.keys()!=expected.keys():
            raise AssertionError(f'{location}: dictionary keys differ')
        for key in expected:
            compare(actual[key],expected[key],f'{location}.{key}')
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(actual)!=len(expected):
            raise AssertionError(f'{location}: list lengths differ')
        for index,(a,e) in enumerate(zip(actual,expected)):
            compare(a,e,f'{location}[{index}]')
    elif actual!=expected:
        raise AssertionError(f'{location}: {actual!r} != {expected!r}')


def verify_qubit_certificates():
    result={'convention':'E00,E01,E10,E11; words multiply left to right',
            'certificates':[qubit_minimum.certificate(m) for m in (1,2,3)]}
    compare(result,load_result('qubit_minimum_results.json'))
    return [{'degree':row['degree'],'rank':row['rank'],'determinant':row['determinant']}
            for row in result['certificates']]


def verify_qutrit_construction():
    result=qutrit_degree8.verify()
    compare(result,load_result('qutrit_results.json'))
    source=load_result('degree8_coefficients.json')
    generated=qutrit_degree8.generate_coefficients_capelli()
    exported={tuple(x-1 for x in term['word']):term['coefficient']
              for term in source['terms']}
    assert generated==exported,'Stored qutrit coefficient dictionary differs'
    assert len(source['terms'])==len(generated),'Duplicate exported words'
    assert source['degree']==8 and source['dimension']==3
    assert source['support_size']==source['l1_norm']==2880
    result['committed_coefficient_export_agrees']=True
    return result


def verify_depolarizing_histograms(exhaustive):
    expected=load_result('depolarizing_results.json')
    qubit=json_value(depolarizing.record(depolarizing.QUBIT_TERMS,2))
    compare(qubit,expected['qubit'])
    assert Fraction(qubit['success_probability_exact'])==Fraction(3,64)

    # The default suite checks the stored exact all-word certificate and also
    # independently recomputes the qutrit probability from the short trace formula.
    source=load_result('degree8_coefficients.json')
    terms=[(tuple(term['word']),term['coefficient']) for term in source['terms']]
    qutrit=expected['qutrit']
    histogram={int(k):v for k,v in qutrit['cycle_histogram'].items()}
    l1=sum(abs(c) for _,c in terms)
    assert qutrit['degree']==8 and qutrit['dimension']==3
    assert qutrit['support_size']==len(terms)==2880
    assert qutrit['l1_norm']==l1==2880
    assert histogram[9]==sum(c*c for _,c in terms)
    assert sum(histogram.values())==sum(c for _,c in terms)**2
    for d in (1,2):
        assert depolarizing.exact_probability(histogram,8,l1,d)==0
    assert depolarizing.exact_probability(histogram,8,l1,3)==Fraction(1,43740)

    if exhaustive:
        assert depolarizing.histogram_cpp(depolarizing.QUBIT_TERMS)==depolarizing.histogram_python(depolarizing.QUBIT_TERMS)
        actual=depolarizing.record(terms,3,use_cpp=True)
        actual['coefficient_source']='degree8_coefficients.json'
        compare(json_value(actual),qutrit)
    return {'qubit_success_exact':qubit['success_probability_exact'],
            'qutrit_success_exact':qutrit['success_probability_exact'],
            'qutrit_cycle_histogram':qutrit['cycle_histogram'],
            'exhaustive_qutrit_pairs_recomputed':exhaustive}


def run_exact_program(program, result_name, output_dir):
    """Run a certificate in isolation, leaving committed records unchanged."""
    destination=output_dir/Path(program).stem
    command=[sys.executable,str(HERE/program),'--output-dir',str(destination)]
    completed=subprocess.run(command,capture_output=True,text=True)
    if completed.returncode:
        raise AssertionError(f'{program} failed:\n{completed.stdout}\n{completed.stderr}')
    return json.loads((destination/result_name).read_text())


def verify_encoding_bound(output_dir):
    prefix='encoding_bound/'
    result=run_exact_program(prefix+'encoding_certificate/verify_certificate.py',
                             'exact_verification.json',output_dir)
    coverage=run_exact_program(prefix+'encoding_orbits/verify_coverage.py',
                               'coverage_status.json',output_dir)
    reduction=run_exact_program(prefix+'encoding_certificate_audit/verify_reduction.py',
                                'reduction_audit.json',output_dir)
    baseline=run_exact_program(prefix+'encoding_certificate_audit/verify_baseline.py',
                               'baseline_audit.json',output_dir)
    assert result['status']==coverage['status']==reduction['status']==baseline['status']=='PASS'
    assert result['certificate_count']==97 and result['exact_positive_matrix_checks']==388
    assert coverage['raw_encodings_covered']==20160 and coverage['combined_classes']==97
    assert result['common_support_bound']=='27/10'
    assert baseline['expected_calls_per_acceptance']=='101/32'
    return {'certificate_count':result['certificate_count'],
            'exact_positive_matrix_checks':result['exact_positive_matrix_checks'],
            'common_support_bound':result['common_support_bound'],
            'minimum_exact_DD_margin':result['minimum_exact_DD_margin'],
            'uniform_operator_perturbation_bound':result['uniform_operator_perturbation_bound'],
            'raw_encodings_covered':coverage['raw_encodings_covered'],
            'cube_classes':coverage['cube_classes'],
            'combined_classes':coverage['combined_classes'],
            'independent_Burnside_count':coverage['independent_Burnside_count'],
            'baseline':baseline,
            'candidate_SHA256':result['candidate_SHA256'],
            'factor_SHA256':result['factor_SHA256']}


def verify_local_controller(output_dir):
    prefix='local_controller/'
    spectral=run_exact_program(prefix+'verify_spectral.py','audit.json',output_dir)
    flags=run_exact_program(prefix+'exact_flags.py','pairwise_exact_flags.json',output_dir)
    thresholds=run_exact_program(prefix+'certify_thresholds.py','certified_thresholds.json',output_dir)
    alternative=run_exact_program(prefix+'physical_proportional_certificate.py',
                                  'physical_proportional_certificate.json',output_dir)
    # Coefficients and the constructive witness must agree exactly with the
    # committed rational records; random checks only supplement the proof.
    compare(flags,load_result(prefix+'pairwise_exact_flags.json'))
    compare(alternative,load_result(prefix+'physical_proportional_certificate.json'))
    compare(thresholds,load_result(prefix+'certified_thresholds.json'))
    assert spectral['antisymmetric_cross_kernel_nonzeros']==12
    assert spectral['antisymmetric_cross_kernel_Hamming_support']==2
    assert alternative['achieved_error']=='3/8'
    assert Fraction(alternative['timecost_in_tau_units'])<Fraction(alternative['baseline12'])
    return {'spectral':spectral,
            'exact_flag_cubic_coefficients_agree':True,
            'thresholds':thresholds,
            'alternative_encoding':alternative}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exhaustive',action='store_true',
                        help='Also compile C++17 and recompute every qutrit word pair')
    parser.add_argument('--report',type=Path,help='Optional machine-readable report path')
    args=parser.parse_args()
    checks={}

    checks['qubit_minimum']=verify_qubit_certificates()
    print('PASS: exact qubit degree-one, degree-two, and degree-three certificates',flush=True)

    checks['qutrit_construction']=verify_qutrit_construction()
    print('PASS: two exact qutrit coefficient generators, source witness, and matrix checks',flush=True)

    switch=qubit_switch_checks.verify()
    compare(json_value(switch),load_result('qubit_switch_results.json'))
    residual_names=('switch_tp_error','identity_choi_error','reference_state_error',
                    'choi_joint_projection_error','input_independence_effect_error',
                    'maximum_scalar_kraus_residual')
    checks['actual_switch']={
        'case_probabilities':{name:row['success_probability'] for name,row in switch.items()
                              if 'success_probability' in row},
        'maximum_residuals':{key:max(row.get(key,0.) for row in switch.values()) for key in residual_names},
        'full_prepared_switch_gauge_error':switch['rectangular_kraus_gauge_change']['complete_prepared_switch_gauge_error'],
        'coefficient_phase_error':switch['coefficient_phase_check']['branch_error'],
    }
    print('PASS: actual qubit switch, nonunital probabilities, references, and Kraus gauge',flush=True)

    access=check_access_model.verify()
    compare(json_value(access),load_result('access_model_checks.json'))
    checks['access_model']=access
    print('PASS: mixed ancillary controls, general effects, classical bound, and correlated noise',flush=True)

    trace=qutrit_trace_crosscheck.verify()
    compare(json_value(trace),load_result('qutrit_trace_results.json'))
    checks['qutrit_trace_probability']=trace
    print('PASS: independent exact 96-term qutrit trace contraction gives 1/43740',flush=True)

    checks['depolarizing_histograms']=verify_depolarizing_histograms(args.exhaustive)
    print('PASS: exact depolarizing histograms'+(' including exhaustive qutrit pairs' if args.exhaustive else ' and stored qutrit certificate'),flush=True)

    regev=regev_probabilities.verify()
    compare(json_value(regev),load_result('regev_results.json'))
    checks['regev_probabilities']=regev
    print('PASS: exact Regev constants, probabilities, and resource figures',flush=True)

    checks['regev_identity']=regev_identity_checks.verify()
    print('PASS: exact Regev matrix-unit evaluations, random tuples, and Cauchy--Binet normalization',flush=True)

    with tempfile.TemporaryDirectory(prefix='quantum-order-certificates-') as temporary:
        output_dir=Path(temporary)
        checks['all_encoding_cost_bound']=verify_encoding_bound(output_dir)
        print('PASS: exact 20,160-encoding coverage, 388 positive-matrix checks, and 101/32 cost witness',flush=True)
        checks['local_controller']=verify_local_controller(output_dir)
        print('PASS: exact pairwise local-noise spectral identities, flag maps, thresholds, and alternative encoding',flush=True)

    report={'passed':True,'python_version':sys.version.split()[0],
            'numpy_version':np.__version__,
            'random_seeds':{'probability_and_construction':20261028,
                            'independent_access_model':202610280314,
                            'local_controller_supplementary':29102026},
            'checks':checks}
    if args.report:
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(f'All {len(checks)} verification groups passed.',flush=True)


if __name__=='__main__':
    main()
