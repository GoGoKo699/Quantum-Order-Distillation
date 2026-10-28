"""Independent exact p_Delta check from Drensky--Kasparian's scalar trace formula.

Uses 96 trace-permutation terms from printed page 752, rather than the 2880
matrix-product words. For dimension 3 the scalar equals the central polynomial.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import argparse
import json


def parity(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))


def trace_terms():
    terms=Counter()
    for p in permutations(range(4)):
        for a,sa in [((4,5),1),((5,4),-1)]:
            for b,sb in [((6,7),1),((7,6),-1)]:
                permutation=[0]*8
                for cycle in [(p[0],p[1])+a,(p[2],p[3])+b]:
                    for x,y in zip(cycle,cycle[1:]+cycle[:1]):
                        permutation[x]=y
                terms[tuple(permutation)]+=2*parity(p)*sa*sb
    return dict(terms)


def verify():
    terms=trace_terms(); histogram=Counter()
    for sigma,cs in terms.items():
        for tau,ct in terms.items():
            inv=[0]*8
            for x,y in enumerate(tau): inv[y]=x
            seen=set(); count=0
            for x in range(8):
                if x in seen: continue
                count+=1; j=x
                while j not in seen:
                    seen.add(j); j=sigma[inv[j]]
            histogram[count]+=cs*ct
    def probability(d):
        return sum((Fraction(weight*d**cycles,d**8*2880**2)
                    for cycles,weight in histogram.items()),Fraction())
    p=probability(3)
    assert len(terms)==96
    assert p==Fraction(1,43740)
    assert probability(1)==probability(2)==0
    result={'source':'Drensky--Kasparian (1985), printed page 752',
            'trace_term_count':len(terms),'trace_cycle_histogram':dict(sorted(histogram.items())),
            'qutrit_success_exact':str(p),'qutrit_success_float':float(p)}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Optional JSON report path')
    args=parser.parse_args()
    text=json.dumps(verify(),indent=2)+'\n'
    print(text,end='')
    if args.output:
        args.output.write_text(text)


if __name__=='__main__': main()
