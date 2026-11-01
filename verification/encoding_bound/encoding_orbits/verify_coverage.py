"""Exact finite coverage certificate. Standard library only.

The integer codewords in rational candidate JSON (preferred) or NPZ solver
output are read without numpy. This file certifies class coverage and
symmetry, not positive semidefiniteness of duals.
"""
from pathlib import Path
from itertools import permutations
import argparse,ast,json,struct,zipfile,hashlib,math
from enumerate_orbits import (ORDERS,SLOT_ACTIONS,COORD_MAPS,cube_canonical,
                              full_canonical,distance,kernel_certificate)

ROOT=Path(__file__).resolve().parent

def cycle_counts(p):
    seen=set();counts={}
    for start in range(len(p)):
        if start in seen:continue
        cur=start;length=0
        while cur not in seen:
            seen.add(cur);length+=1;cur=p[cur]
        counts[length]=counts.get(length,0)+1
    return counts

def burnside_count():
    cube_actions=[tuple(mapping[v^shift] for v in range(8)) for mapping in COORD_MAPS for shift in range(8)]
    assert len(set(cube_actions))==48
    totals={1:0,2:0,3:0};contributions=[]
    for p in SLOT_ACTIONS:
        dom=cycle_counts(p);order=max(dom)
        for a in cube_actions:
            ran=cycle_counts(a);fixed=1
            for length,n in dom.items():
                count=ran.get(length,0)
                if count<n:fixed=0;break
                fixed*=math.factorial(count)//math.factorial(count-n)*length**n
            totals[order]+=fixed
            if fixed:contributions.append({'domain_cycle_length':order,'cube_cycles':ran,'fixed_encodings':fixed})
    assert totals=={1:20160,2:7488,3:288}
    total=sum(totals.values());assert total==27936 and total//288==97
    return {'identity_domain_contribution':totals[1],'transposition_domain_contribution':totals[2],
            'three_cycle_domain_contribution':totals[3],'fixed_point_sum':total,'group_order':288,'orbit_count':97,
            'free_cube_involutions':sum(cycle_counts(a)=={2:4} for a in cube_actions),
            'cube_order_three_actions':sum(cycle_counts(a)=={1:2,3:2} for a in cube_actions)}

def read_npz_codes(path):
    with zipfile.ZipFile(path) as z:data=z.read('codes.npy')
    assert data[:6]==b'\x93NUMPY'
    major,minor=data[6:8]
    if major==1:n=struct.unpack('<H',data[8:10])[0];offset=10
    elif major in (2,3):n=struct.unpack('<I',data[8:12])[0];offset=12
    else:raise AssertionError('Unsupported NPY format')
    header=ast.literal_eval(data[offset:offset+n].decode('latin1'))
    assert header['shape']==(6,) and not header['fortran_order']
    descr=header['descr'];assert descr in ('<i8','<i4','<i2','|i1','<u8','<u4','<u2','|u1')
    widths={'8':8,'4':4,'2':2,'1':1};width=widths[descr[-1]]
    raw=data[offset+n:];assert len(raw)==6*width
    codes=tuple(int.from_bytes(raw[i*width:(i+1)*width],'little',signed='i' in descr) for i in range(6))
    assert len(set(codes))==6 and all(0<=v<8 for v in codes)
    return codes

def witness(codes):
    best=None
    for si,p in enumerate(SLOT_ACTIONS):
        permuted=tuple(codes[i] for i in p)
        shift=permuted[0]
        for ci,mapping in enumerate(COORD_MAPS):
            transformed=tuple(mapping[c^shift] for c in permuted)
            row=(transformed,si,shift,ci)
            if best is None or row<best:best=row
    assert best[0]==full_canonical(codes)
    return best

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--certificate-dir',default=str(ROOT.parent/'encoding_dual'))
    ap.add_argument('--certificate-json',default=None)
    ap.add_argument('--output-dir',type=Path,default=ROOT)
    args=ap.parse_args();certificate_dir=Path(args.certificate_dir)
    output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=True)
    cubes=sorted({cube_canonical((0,)+p) for p in permutations(range(1,8),5)})
    classes=sorted({full_canonical(c) for c in cubes})
    assert len(cubes)==420 and len(classes)==97
    indices={c:i for i,c in enumerate(classes)}
    compact=[];hist=[0]*97
    for c in cubes:
        rep,si,shift,ci=witness(c)
        compact.append({'cube_representative':list(c),'combined_class':indices[rep],
                        'slot_action':list(SLOT_ACTIONS[si]),'xor_translation':shift,
                        'coordinate_map':list(COORD_MAPS[ci])})
    for c in permutations(range(8),6):
        rep,si,shift,ci=witness(c)
        reconstructed=tuple(COORD_MAPS[ci][c[j]^shift] for j in SLOT_ACTIONS[si])
        assert reconstructed==rep
        hist[indices[rep]]+=1
    assert sum(hist)==20160
    candidate=Path(args.certificate_json) if args.certificate_json else certificate_dir/'certificate_candidates.json'
    certs=[]
    if candidate.is_file():
        data=json.loads(candidate.read_text());rows=data['certificates']
        assert len(rows)==97
        digest=hashlib.sha256(candidate.read_bytes()).hexdigest()
        for i,row in enumerate(rows):
            codes=tuple(row['codes'])
            assert len(codes)==6 and all(type(c) is int and 0<=c<8 for c in codes) and len(set(codes))==6
            certs.append((candidate.name+'#'+str(i),codes,digest))
        source_format='rational_candidate_json'
    else:
        paths=sorted(certificate_dir.glob('dual_*.npz'))
        certs=[(p.name,read_npz_codes(p),hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]
        source_format='npz_integer_codes'
    assert len(certs)==97, f'Expected 97 dual certificates, got {len(certs)}'
    coverage={};records=[]
    for source,codes,digest in certs:
        rep=full_canonical(codes);idx=indices[rep]
        assert idx not in coverage, f'Duplicate class {idx}'
        coverage[idx]=source
        records.append({'certificate_source':source,'codes':list(codes),'class_index':idx,
                        'canonical_representative':list(rep),
                        'sha256':digest})
    assert set(coverage)==set(range(97))
    kernel=kernel_certificate()
    (output/'cube_to_combined_map.json').write_text(json.dumps(compact,separators=(',',':'))+'\n')
    result={'status':'PASS','raw_encodings_covered':20160,'cube_classes':420,
            'combined_classes':97,'saved_certificates':len(certs),'certificate_format':source_format,'duplicate_or_missing_classes':0,
            'orbit_size_histogram':{str(n):hist.count(n) for n in sorted(set(hist))},
            'exact_transform_reconstructions':20160,'kernel_certificate':kernel,'independent_Burnside_count':burnside_count(),
            'scope':'Exact encoding coverage and symmetry only; numerical/rational dual validity is checked separately.',
            'dual_mapping':records}
    (output/'coverage_status.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='dual_mapping'},indent=2))

if __name__=='__main__':main()
