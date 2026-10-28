"""Exact normalization and asymptotic resource figures for the Regev protocol.

Central-polynomial identity: Procesi (2021), A note on the Formanek Weingarten
function, Eq.(1.65), Theorem 7, Eq.(2.85), pp.94,99--100,
DOI 10.1285/i15900932v41n1p69. Also in arXiv:2102.10657v2, Theorem 4.13.
This script evaluates the inherited constants.
"""
from fractions import Fraction
from math import factorial, lgamma, log
from pathlib import Path
import argparse
import json


def constant(d):
    c=Fraction(1)
    for k in range(1,d+1): c*=factorial(2*k-1)
    for k in range(1,d): c/=factorial(k)
    a=(-1)**(d-1)*c*c/(factorial(d)**2*(2*d-1))
    return c,a


def log_probability(d):
    logc=sum(lgamma(2*k) for k in range(1,d+1))-sum(lgamma(k+1) for k in range(1,d))
    return (4*logc-4*lgamma(d+1)-2*log(2*d-1)
            -2*lgamma(d*d+1)-2*d*d*log(d))


def verify():
    rows=[]
    for d in (2,3,4,5):
        c,a=constant(d); r=d*d; support=factorial(r)**2
        p=a*a/(factorial(r)**2*d**(2*r))
        row={'dimension':d,'queries':2*r,'support_size':support,
             'control_qubits':(support-1).bit_length(),'C_d':str(c),'a_d':str(a),
             'depolarizing_success_exact':str(p),'depolarizing_success_float':float(p),
             'ln_success':log_probability(d),'expected_query_count_exact':str(Fraction(2*r,1)/p)}
        rows.append(row)
    assert rows[0]['a_d']=='-3' and rows[0]['depolarizing_success_exact']=='1/16384'
    assert rows[1]['a_d']=='720' and rows[1]['depolarizing_success_exact']=='1/98411002933824'
    asymptotic=[]
    for d in (10,30,100,300):
        leading=-4*d*d*log(d)+(4*log(2)-1)*d*d
        refined=leading-2*d*log(d)+(2*log(2)+2)*d-(35/6)*log(d)
        asymptotic.append({'dimension':d,'ln_success':log_probability(d),
                           'leading_residual_over_d_log_d':(log_probability(d)-leading)/(d*log(d)),
                           'refined_residual':log_probability(d)-refined})
    return {'exact_examples':rows,'asymptotic_checks':asymptotic}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Optional JSON report path')
    args=parser.parse_args()
    text=json.dumps(verify(),indent=2)+'\n'
    print(text,end='')
    if args.output:
        args.output.write_text(text)


if __name__=='__main__': main()
