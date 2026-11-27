"""Independent unpruned enumeration, using polynomial arithmetic (no log tables)."""
from itertools import product
from math import gcd
from pathlib import Path
import json, subprocess

def field(N,p):
    def mul(a,b):
        z=0
        while b:
            if b&1:z^=a
            b>>=1;a<<=1
            if a&(1<<N):a^=p
        return z
    return mul

def rank(V):
    B={}
    for x in V:
        while x:
            k=x.bit_length()-1
            if k in B:x^=B[k]
            else:B[k]=x;break
    return len(B)

def spaces(N,n):
    def rec(B,piv):
        if len(B)==n:yield B;return
        for p in range(piv[-1]+1,N-(n-len(B))+1):
            free=[i for i in range(1,p) if i not in piv]
            for mask in range(1<<len(free)):
                x=(1<<p)|sum((1<<b) for k,b in enumerate(free) if mask>>k&1)
                yield from rec(B+[x],piv+[p])
    yield from rec([1],[0])

def span(B):
    W={0}
    for x in B:W|={w^x for w in W}
    return W

def rowmul(A,B):
    return tuple(__import__('functools').reduce(int.__xor__,
        (b for i,b in enumerate(B) if a>>i&1),0) for a in A)

def coordinates(B,x):
    E={}
    for j,a in enumerate(B):
        z=1<<j
        while a:
            k=a.bit_length()-1
            if k in E:a^=E[k][0];z^=E[k][1]
            else:E[k]=(a,z);break
    z=0
    while x:
        k=x.bit_length()-1;x^=E[k][0];z^=E[k][1]
    return z

def feedback(B,n,m,mul):
    V=[];pw=1
    for j in range(m):
        V += [mul(pw,b) for b in B];pw=mul(pw,2)
    c=[coordinates(V,mul(pw,b)) for b in B]
    return [tuple((x>>(j*n))&((1<<n)-1) for x in c) for j in range(m)]

def verify(n,m,p,centralizer=True):
    N=n*m;Q=1<<N;mul=field(N,p)
    powers=[];x=1
    for j in range(Q-1):powers.append(x);x=mul(x,2)
    assert x==1 and len(set(powers))==Q-1
    primitive=[powers[j] for j in range(Q-1) if gcd(j,Q-1)==1]
    def independent(B,beta):
        V=[];pw=1
        for j in range(m):
            V += [mul(pw,b) for b in B];pw=mul(pw,beta)
        return rank(V)==N
    total=good=split=0;good_sets=[]
    for B in spaces(N,n):
        total+=1;W=span(B)
        if all(independent(B,beta) for beta in primitive):
            good+=1;good_sets.append(sorted(W))
        if not independent(B,2):continue
        split+=1
        if not centralizer:continue
        A=feedback(B,n,m,mul)
        commuting=all(rowmul(a,b)==rowmul(b,a) for a in A for b in A)
        L=[c for c in range(Q) if all(mul(c,b) in W for b in B)]
        C=[]
        for z in range(1<<(n*n)):
            mat=tuple((z>>(i*n))&((1<<n)-1) for i in range(n))
            if all(rowmul(mat,a)==rowmul(a,mat) for a in A):C.append(mat)
        assert len(C)==len(L)
        assert commuting==(len(L)==(1<<n))
        # Independently form the subspace polynomial as product (X+w).
        f=[1]
        for w in W:
            g=[0]*(len(f)+1)
            for i,c in enumerate(f):g[i]^=mul(c,w);g[i+1]^=c
            f=g
        assert all(c==0 or (i and i&(i-1)==0) for i,c in enumerate(f))
        r=n
        for i in range(n):
            if f[1<<i]:r=gcd(r,n-i)
        assert len(L)==1<<r
    E=[]
    for x in range(Q):
        y=x
        for _ in range(n):y=mul(y,y)
        if y==x:E.append(x)
    assert good_sets==[E]
    return dict(n=n,m=m,polynomial=hex(p),unpruned_spaces=total,
                alpha_splitting=split,universal_spaces=good,
                centralizer_and_subspace_polynomial_checked=centralizer)

if __name__=='__main__':
    results=[verify(2,2,0x13),verify(2,3,0x43),verify(3,2,0x43),verify(4,2,0x11d,False)]
    base=Path(__file__).parent
    for r in results:
        exe=base/('search_quadratic' if r['m']==2 else 'search_general')
        args=[str(exe),str(r['n'])]
        if r['m']!=2:args.append(str(r['m']))
        args.append(r['polynomial'][2:])
        out=subprocess.check_output(args,text=True)
        assert out.count('GOOD ')==r['universal_spaces']
        assert 'is_E=1' in out
    (base/'independent_verification.json').write_text(json.dumps(results,indent=2))
    print(json.dumps(results,indent=2))
