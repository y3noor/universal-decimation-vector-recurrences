"""Verify the finite feedback enumerations."""
from itertools import product
from math import gcd
import json
from pathlib import Path

def rank(rows):
    piv = {}
    for x in rows:
        while x:
            b=x.bit_length()-1
            if b in piv: x ^= piv[b]
            else:
                piv[b]=x
                break
    return len(piv)

def mul(A,B):
    out=[]
    for a in A:
        r=0
        for i,b in enumerate(B):
            if (a>>i)&1:r ^= b
        out.append(r)
    return tuple(out)

def identity(N):return tuple(1<<i for i in range(N))
def power(A,d):
    b=identity(len(A))
    while d:
        if d&1:b=mul(b,A)
        A=mul(A,A);d>>=1
    return b

def apply(A,x):
    return sum(((a&x).bit_count()&1)<<i for i,a in enumerate(A))

def companion(blocks,n):
    m=len(blocks);N=n*m
    rows=[1<<(i+n) for i in range(N-n)]
    rows += [sum(blocks[j][r]<<(j*n) for j in range(m)) for r in range(n)]
    return tuple(rows)

def full_cycle(A):
    N=len(A);tau=(1<<N)-1
    x=1;seen=set()
    for _ in range(tau):
        if x==0 or x in seen:return False
        seen.add(x);x=apply(A,x)
    return x==1 and len(seen)==tau

def obs_rank(A,n,m,d):
    P=identity(n); P=tuple(P)  # projection rows in N-bit representation
    T=power(A,d);rows=[]
    for _ in range(m):
        rows.extend(P);P=mul(P,T)
    return rank(rows)

def commute(blocks):
    return all(mul(A,B)==mul(B,A) for A in blocks for B in blocks)

def classify(n,m):
    mats=[tuple((x>>(r*n))&((1<<n)-1) for r in range(n))
          for x in range(1<<(n*n))]
    N=n*m;tau=(1<<N)-1;ds=[d for d in range(1,tau) if gcd(d,tau)==1]
    result=dict(n=n,m=m,examined=0,maximal=0,commuting=0,
                admissible=0,admissible_noncommuting=0,coprime_decimations=len(ds))
    for blocks in product(mats,repeat=m):
        result['examined']+=1
        A=companion(blocks,n)
        if not full_cycle(A):continue
        result['maximal']+=1
        c=commute(blocks);result['commuting']+=c
        ok=all(obs_rank(A,n,m,d)==N for d in ds)
        result['admissible']+=ok
        result['admissible_noncommuting']+=ok and not c
    return result

if __name__=='__main__':
    results=[classify(2,2),classify(2,3)]
    assert results[0]['examined']==256 and results[0]['maximal']==16
    assert results[0]['admissible']==4 and results[0]['admissible_noncommuting']==0
    assert results[1]['examined']==4096 and results[1]['maximal']==192
    assert results[1]['admissible']==12 and results[1]['admissible_noncommuting']==0
    blocks=((2,1),(0,0),(1,0))
    A=companion(blocks,2)
    assert full_cycle(A) and not commute(blocks)
    assert obs_rank(A,2,3,2)==6
    failed=[d for d in range(1,63) if gcd(d,63)==1 and obs_rank(A,2,3,d)<6]
    assert failed
    report=dict(finite_classifications=results,
                obstruction_counterexample={'blocks':blocks,'state_period':63,
                 'rank_O2':6,'failing_decimations':failed},
                general_question='not proved')
    Path(__file__).with_name('results.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
