"""Exact computational-basis encoding quotient for three slots, stdlib only."""
from itertools import permutations, product
from fractions import Fraction as F
from pathlib import Path
import json, hashlib
ROOT=Path(__file__).resolve().parent
ORDERS=tuple(permutations(range(3)))
COORDS=tuple(permutations(range(3)))
SLOT_ACTIONS=tuple(tuple(ORDERS.index(tuple(g[x] for x in o)) for o in ORDERS) for g in permutations(range(3)))

def coordinate_image(code, perm):
    return sum(((code>>perm[j])&1)<<j for j in range(3))
COORD_MAPS=tuple(tuple(coordinate_image(c,p) for c in range(8)) for p in COORDS)

def cube_canonical(codes):
    # Cube isometry = XOR translation followed by coordinate permutation.
    # Sending c0 to zero removes the translation freedom uniquely.
    x=tuple(c^codes[0] for c in codes)
    return min(tuple(mapping[c] for c in x) for mapping in COORD_MAPS)

def full_canonical(codes):
    return min(cube_canonical(tuple(codes[i] for i in p)) for p in SLOT_ACTIONS)

def distance(codes):
    return tuple(tuple((a^b).bit_count() for b in codes) for a in codes)

def mm(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]

def kernel_certificate():
    mats=[[[1,0],[1,0]],[[0,1],[0,-1]]]
    histories=[]
    for hist in product(range(2),repeat=3):
        rows=[]
        for o in ORDERS:
            a=[[1,0],[0,1]]
            for k in o:a=mm(a,mats[hist[k]])
            rows.append(tuple(v for row in a for v in row))
        histories.append(rows)
    # Kraus products share 1/sqrt(8), normalized Choi vectors 1/sqrt(2).
    # Rational target coordinate kernels = integer Gram /16.
    ks=[]
    for a,b in product(range(4),repeat=2):
        k=tuple(tuple(F(sum(h[i][a]*h[j][b] for h in histories),16) for j in range(6)) for i in range(6))
        for p in SLOT_ACTIONS:
            assert all(k[i][j]==k[p[i]][p[j]] for i,j in product(range(6),repeat=2))
        ks.append(k)
    # Individual matrix entries include the U--V target cross kernel.
    # Hence equality is stronger than a check of the four Bell diagonals.
    return {'target_coordinate_kernels':len(ks),'slot_actions':len(SLOT_ACTIONS),
            'exact_entry_invariances':len(ks)*len(SLOT_ACTIONS)*36,
            'result':'all 16 full target-Choi coordinate kernels invariant under left physical-slot relabeling'}

def main():
    cubereps=sorted({cube_canonical((0,)+p) for p in permutations(range(1,8),5)})
    assert len(cubereps)==420
    representatives=sorted({full_canonical(c) for c in cubereps})
    repindex={c:i for i,c in enumerate(representatives)}
    rows=[];hist=[0]*len(representatives);cubehist=[0]*len(representatives)
    for c in cubereps:cubehist[repindex[full_canonical(c)]]+=1
    for c in permutations(range(8),6):
        rep=full_canonical(c);r=repindex[rep];hist[r]+=1
        rows.append({'codewords':list(c),'representative_index':r})
    assert len(rows)==20160 and sum(hist)==20160
    assert all(a==48*b for a,b in zip(hist,cubehist))
    # Direct metric equivalence certificate via the generating group actions.
    checks=0
    for c in cubereps:
        d=distance(c)
        for mapping in COORD_MAPS:
            transformed=tuple(mapping[v] for v in c)
            assert distance(transformed)==d;checks+=1
        for shift in range(8):
            assert distance(tuple(v^shift for v in c))==d;checks+=1
        for p in SLOT_ACTIONS:
            dp=distance(tuple(c[i] for i in p))
            assert all(dp[i][j]==d[p[i]][p[j]] for i,j in product(range(6),repeat=2));checks+=1
    data=[]
    for c,n,cn in zip(representatives,hist,cubehist):
        unused=sorted(set(range(8))-set(c))
        data.append({'index':repindex[c],'codes':list(c),'codewords':[format(v,'03b') for v in c],
                     'distance_matrix':[list(r) for r in distance(c)],'full_orbit_size':n,
                     'cube_orbits_in_class':cn,'unused_vertex_distance':(unused[0]^unused[1]).bit_count(),
                     'joint_stabilizer_size':288//n})
    out={'status':'PASS','raw_injective_encodings':20160,'cube_group_size':48,
         'cube_orbits':len(cubereps),'slot_group_size':6,'combined_group_size':288,
         'combined_orbits':len(representatives),'orbit_size_histogram':{str(n):hist.count(n) for n in sorted(set(hist))},
         'cube_metric_generating_checks':checks,'kernel_certificate':kernel_certificate(),'representatives':data}
    (ROOT/'representatives.json').write_text(json.dumps(out,indent=2)+'\n')
    (ROOT/'all_encoding_orbit_map.json').write_text(json.dumps(rows,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='representatives'},indent=2))
    print('Map SHA256',hashlib.sha256((ROOT/'all_encoding_orbit_map.json').read_bytes()).hexdigest())
if __name__=='__main__':main()
