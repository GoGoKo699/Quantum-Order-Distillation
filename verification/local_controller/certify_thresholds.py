from fractions import Fraction as F
from pathlib import Path
import json,math
import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent,
                     help='Directory for regenerated exact certificate records')
_args=_parser.parse_args() if __name__=='__main__' else _parser.parse_args([])
ROOT=_args.output_dir.resolve()
ROOT.mkdir(parents=True,exist_ok=True)

def ln09_bounds(n=45):
 z=F(1,19);s=-2*sum(z**(2*k+1)/F(2*k+1) for k in range(n));tail=2*z**(2*n+1)/((2*n+1)*(1-z*z));return s-tail,s

def expneg_bounds(v,n=80):
 assert -1<v<0
 # n even => upper; n+1 odd =>lower
 s=sum(v**k/F(math.factorial(k)) for k in range(n+1));return s+v**(n+1)/math.factorial(n+1),s

def q_bounds(r):
 lo,hi=ln09_bounds();return expneg_bounds(r*lo)[0],expneg_bounds(r*hi)[1]

def h(q):return (-106+69*q+69*q*q+92*q**3)/64

def pquery(q):return 92*q**3+69*q*q+69*q-202

def proportional(t):return 4*t**11+24*t**9+3*t**8+18*t**6+3*t**5+18*t**3-54*t*t-12

def bracket_decimal(fun,lo,hi):
 a=F(lo);b=F(hi);va=fun(a);vb=fun(b);assert va<0<vb;return {'lower':str(a),'upper':str(b),'lower_decimal':lo,'upper_decimal':hi,'sign_lower':-1,'sign_upper':1,'lower_value_float':float(va),'upper_value_float':float(vb)}

if __name__=='__main__':
 rlo=F('1.14647049186791');rhi=F('1.14647049186792');lolo,lohi=q_bounds(rlo);hilo,hihi=q_bounds(rhi)
 vlo=h(lolo)-rlo;vhi=h(hihi)-rhi;assert vlo>0>vhi
 # This outward rational bracket is checked exactly; no numeric solver is trusted.
 ta=F('0.98809436711299');tb=F('0.98809436711300')
 assert proportional(ta)<0<proportional(tb)
 rt=float((ta+tb)/2);print('troot',rt,'xroot',rt**2,'qroot',rt**3)
 # Exact positive derivative bound for t∈[9/10,1].
 a=F(9,10);deriv_lower=44*a**10+216*a**8+24*a**7+108*a**5+15*a**4+54*a*a-108
 assert deriv_lower>0
 out={'time_ratio_root':{'lo':str(rlo),'hi':str(rhi),'F_lo_lower_bound':float(vlo),'F_hi_upper_bound':float(vhi),'rigor':'Fraction atanhlike series bounds forlog(.9), alternatingexp; h(q)increasing andh(.9^r)-r strictlydecreasing'},'query_q_root':bracket_decimal(pquery,'0.939452711560243','0.939452711560244'),'proportional_root':{'tlo':str(ta),'thi':str(tb),'xlo':str(ta*ta),'xhi':str(tb*tb),'t_float':rt,'x_float':rt*rt,'q_float':rt**3,'gamma_tau':-math.log(rt),'polynomial_coefficients_ascending':[-12,0,-54,18,0,3,18,0,3,24,0,4],'strict_derivative_lower_bound_on_t_09_1':str(deriv_lower)}}
 print(json.dumps(out,indent=2));(ROOT/'certified_thresholds.json').write_text(json.dumps(out,indent=2))
