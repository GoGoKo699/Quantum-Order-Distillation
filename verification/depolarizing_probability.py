"""Exact depolarizing success probability from an order-polynomial coefficient list.

No sampling over Kraus histories: matrix-unit contractions give the exact Gram
entry d**(cycles(C_pi C_tau^{-1})-m-1). The polynomial must be central to interpret
the maximally-mixed-input success probability as input-independent identity success.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import json
import subprocess
import tempfile


QUBIT_TERMS = [
    ((1,2,3,4),1), ((1,2,4,3),-1),
    ((2,1,3,4),-1), ((2,1,4,3),1),
    ((3,4,1,2),1), ((3,4,2,1),-1),
    ((4,3,1,2),-1), ((4,3,2,1),1),
]


def histogram_python(terms):
    m=len(terms[0][0])
    cycles=[]
    inverses=[]
    for word,_ in terms:
        c=[0]*(m+1)
        for a,b in zip((0,)+tuple(word),tuple(word)+(0,)):
            c[a]=b
        inv=[0]*(m+1)
        for a,b in enumerate(c): inv[b]=a
        cycles.append(c); inverses.append(inv)
    histogram={k:0 for k in range(1,m+2)}
    for i,(_,ci) in enumerate(terms):
        for j in range(i,len(terms)):
            cj=terms[j][1]
            visited=set(); count=0
            for k in range(m+1):
                if k in visited: continue
                count+=1; pos=k
                while pos not in visited:
                    visited.add(pos)
                    pos=cycles[i][inverses[j][pos]]
            histogram[count] += (1 if i==j else 2)*ci*cj
    return histogram


def histogram_cpp(terms):
    source=Path(__file__).with_name('cycle_histogram.cpp')
    m=len(terms[0][0])
    payload=f'{len(terms)} {m}\n'+''.join(
        f'{c} '+ ' '.join(map(str,w))+'\n' for w,c in terms)
    with tempfile.TemporaryDirectory(prefix='order_cycles_') as folder:
        executable=Path(folder)/'cycle_histogram'
        subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(executable)],check=True)
        result=subprocess.run([str(executable)],input=payload,text=True,capture_output=True,check=True)
    return {int(k):int(v) for k,v in (line.split() for line in result.stdout.splitlines())}


def exact_probability(histogram, degree, l1, dimension):
    return sum((Fraction(weight*dimension**cycles,dimension**(degree+1)*l1*l1)
                for cycles,weight in histogram.items()),Fraction())


def record(terms,dimension,use_cpp=False):
    m=len(terms[0][0]); l1=sum(abs(c) for _,c in terms)
    histogram=histogram_cpp(terms) if use_cpp else histogram_python(terms)
    assert sum(histogram.values())==sum(c for _,c in terms)**2
    assert histogram[m+1]==sum(c*c for _,c in terms)
    p=exact_probability(histogram,m,l1,dimension)
    return {'degree':m,'dimension':dimension,'support_size':len(terms),'l1_norm':l1,
            'cycle_histogram':{k:v for k,v in histogram.items() if v},
            'success_probability_exact':str(p),'success_probability_float':float(p),
            'success_at_dimension_1':str(exact_probability(histogram,m,l1,1)),
            'success_at_dimension_2':str(exact_probability(histogram,m,l1,2))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coefficients',type=Path,
                        default=Path(__file__).with_name('degree8_coefficients.json'))
    parser.add_argument('--exhaustive',action='store_true',
                        help='Recompute all qutrit word pairs using the C++17 helper')
    parser.add_argument('--output',type=Path,help='Optional JSON report path')
    args=parser.parse_args()
    qubit=record(QUBIT_TERMS,2)
    assert qubit['cycle_histogram']=={1:16,3:-24,5:8}
    assert qubit['success_probability_exact']=='3/64'
    results={'qubit':qubit}
    if args.exhaustive:
        assert histogram_python(QUBIT_TERMS)==histogram_cpp(QUBIT_TERMS)
        data=json.loads(args.coefficients.read_text())
        terms=[(tuple(item['word']),item['coefficient']) for item in data['terms']]
        result=record(terms,data['dimension'],use_cpp=True)
        result['coefficient_source']=args.coefficients.name
        results['qutrit']=result
        assert result['success_at_dimension_1']=='0'
        assert result['success_at_dimension_2']=='0'
        assert result['success_probability_float']>0
    text=json.dumps(results,indent=2)+'\n'
    print(text,end='')
    if args.output: args.output.write_text(text)


if __name__=='__main__': main()
