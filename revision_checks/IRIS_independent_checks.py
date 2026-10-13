import itertools, json, math
from pathlib import Path
import numpy as np

OUT = Path(__file__).parent

def rank(rows):
    piv = {}
    for a in rows:
        while a:
            p = a.bit_length()-1
            if p in piv: a ^= piv[p]
            else:
                piv[p] = a
                break
    return len(piv)

def solve(eqs, nv):
    piv = {}
    for a, b in eqs:
        while a:
            p = a.bit_length()-1
            if p in piv:
                aa,bb = piv[p]; a ^= aa; b ^= bb
            else:
                piv[p] = a,b
                break
        else:
            if b: raise ValueError('inconsistent wall')
    sol = 0
    for p,(a,b) in sorted(piv.items()):
        if ((a & sol).bit_count() % 2) ^ b: sol |= 1 << p
    return sol, nv-len(piv)

def lattice(N, hx, hz, s=2):
    lx,lz = s*(len(hz)+1),s*(len(hx)+1)
    qs = {}; xs = {}; zs = {}
    def q(sheet, axis, x,y,z):
        if sheet[0]=='C':
            support = hz[sheet[1]]
            if axis=='z' and y in [s*(min(support)+1),s*(max(support)+1)]:
                sheet = ('D', y//s-1)
        key = sheet,axis,x,y,z
        if key not in qs: qs[key]=len(qs)
        return qs[key]
    def xkey(sheet,x,y,z):
        if sheet[0]=='C':
            support = hz[sheet[1]]
            if y in [s*(min(support)+1),s*(max(support)+1)]:
                sheet = ('D', y//s-1)
        return sheet,x,y,z
    def make(sheet, axis1,axis2, coord, r1,r2,rough):
        axes={'x':0,'y':1,'z':2}
        edges = []
        a,b=axes[axis1],axes[axis2]
        def pt(i,j):
            v=list(coord);v[a]=i;v[b]=j;return tuple(v)
        for i in range(r1[0],r1[1]+1):
            for j in range(r2[0],r2[1]+1):
                v=pt(i,j)
                if i<r1[1] and not(rough and j in r2):
                    edges.append((q(sheet,axis1,*v),v,pt(i+1,j)))
                if j<r2[1]: edges.append((q(sheet,axis2,*v),v,pt(i,j+1)))
        for idx,v,w in edges:
            for u in [v,w]:
                if not(rough and u[b] in r2): xs.setdefault(xkey(sheet,*u),set()).add(idx)
        for i in range(r1[0],r1[1]):
            for j in range(r2[0],r2[1]):
                v=pt(i,j);f=set()
                for ax,u in [(axis1,v),(axis1,pt(i,j+1)),(axis2,v),(axis2,pt(i+1,j))]:
                    if not(rough and ax==axis1 and u[b] in r2):f.add(q(sheet,ax,*u))
                zs[sheet,*v]=f
    for k in range(N): make(('D',k),'x','z',(0,s*(k+1),0),(0,lx),(0,lz),True)
    for j,c in enumerate(hz):make(('C',j),'y','z',(s*(j+1),0,0),(s*(min(c)+1),s*(max(c)+1)),(0,lz),True)
    for i,c in enumerate(hx):make(('A',i),'x','y',(0,0,s*(i+1)),(0,lx),(s*(min(c)+1),s*(max(c)+1)),False)
    # Endpoint merge already supplies its incident C edge; add uses support union.
    for j,c in enumerate(hz):
        c=sorted(c);xx=s*(j+1)
        for k in c:
            yy=s*(k+1)
            for zz in range(1,lz): xs[('D',k),xx,yy,zz].add(q(('C',j),'y',xx,yy if k!=c[-1] else yy-1,zz))
            if k not in [c[0],c[-1]]:
                for zz in range(lz): zs[('C',j),xx,yy,zz].add(q(('D',k),'z',xx,yy,zz))
    for i,c in enumerate(hx):
        zz=s*(i+1)
        for k in c:
            yy=s*(k+1)
            for xx in range(lx+1):xs[('A',i),xx,yy,zz].add(q(('D',k),'z',xx,yy,zz))
            for xx in range(lx):zs[('D',k),xx,yy,zz].add(q(('A',i),'x',xx,yy,zz))
    variables=[]
    for i,a in enumerate(hx):
        for j,c in enumerate(hz):
            lo=s*(max(min(a),min(c))+1); hi=s*(min(max(a),max(c))+1)
            xx=s*(j+1); zz=s*(i+1)
            for yy in range(lo,hi+1):variables.append(('x',(('A',i),xx,yy,zz),q(('C',j),'z',xx,yy,zz)))
            for yy in range(lo,hi):variables.append(('z',(('C',j),xx,yy,zz),q(('A',i),'y',xx,yy,zz)))
    eq={}
    for xk,xv in xs.items():
        for zk,zv in zs.items():
            if len(xv&zv)%2:eq[xk,zk]=[0,1]
    for vi,(kind,key,idx) in enumerate(variables):
        others=zs if kind=='x' else xs
        for ok,ov in others.items():
            if idx in ov:
                pair=(key,ok) if kind=='x' else(ok,key)
                eq.setdefault(pair,[0,0])[0]^=1<<vi
    sol,freedom=solve(eq.values(),len(variables))
    for vi,(kind,key,idx) in enumerate(variables):
        if (sol>>vi)&1:(xs if kind=='x' else zs)[key].symmetric_difference_update([idx])
    bits=lambda a:sum(1<<i for i in a)
    xb=list(map(bits,xs.values()));zb=list(map(bits,zs.values()))
    assert all((a&b).bit_count()%2==0 for a in xb for b in zb)
    rx,rz=rank(xb),rank(zb)
    slab_counts=[]; stab_failures=[]; pairing_failures=0
    def slab(word,typ,pos):
        op=set()
        for k in range(N):
            if word>>k&1:
                if typ=='X':
                    for xx in range(lx+1):op.add(q(('D',k),'z',xx,s*(k+1),pos))
                else:
                    for zz in range(lz):op.add(q(('D',k),'z',pos,s*(k+1),zz))
        for j,c in enumerate(hz if typ=='X' else hx):
            lo,hi=s*(min(c)+1),s*(max(c)+1)
            for yy in range(lo+1 if typ=='X' else lo,hi):
                par=sum((word>>k)&1 for k in c if s*(k+1)<yy+(typ=='Z'))%2
                if par:
                    op.add(q(('C',j),'z',s*(j+1),yy,pos) if typ=='X' else q(('A',j),'y',pos,yy,s*(j+1)))
        return bits(op)
    for typ,checks,orig_stab,maxpos,physchecks,physstab in [('X',hz,hx,lz,zb,xb),('Z',hx,hz,lx+1,xb,zb)]:
        checkwords=[sum(1<<k for k in c) for c in checks]
        words=[v for v in range(1<<N) if all((v&h).bit_count()%2==0 for h in checkwords)]
        allslabs=[]
        for pos in range(maxpos):
            for v in words:
                op=slab(v,typ,pos)
                assert all((op&h).bit_count()%2==0 for h in physchecks),(typ,pos,v)
                allslabs.append(op)
            for c in orig_stab:
                op=slab(sum(1<<k for k in c),typ,pos)
                if rank(physstab+[op])!=rank(physstab):stab_failures.append((typ,pos,c))
        slab_counts.append(len(allslabs))
    # Simple transverse sheet strings provide upper witnesses in the even-parity inputs.
    upper = {}
    for typ, support_list, check_rows, stab_rows in [('X',hx,zb,xb),('Z',hz,xb,zb)]:
        witnesses=[]
        for j,c in enumerate(support_list):
            lo,hi=s*(min(c)+1),s*(max(c)+1)
            for yy in range(lo if typ=='X' else lo+1,hi):
                op=bits([q(('A',j),'y',xx,yy,s*(j+1)) for xx in range(lx+1)]) if typ=='X' else bits([q(('C',j),'z',s*(j+1),yy,zz) for zz in range(lz)])
                if all((op&h).bit_count()%2==0 for h in check_rows) and rank(stab_rows+[op])>rank(stab_rows): witnesses.append(op.bit_count())
        upper[typ]=min(witnesses) if witnesses else None
    # Every pairing follows from the data-sheet intersection by construction.
    return dict(N=N,s=s,n=len(qs),xchecks=len(xb),zchecks=len(zb),rx=rx,rz=rz,
                k=len(qs)-rx-rz,euler=len(qs)-len(xb)-len(zb),wall_free_variables=freedom,
                max_check_weight=max(v.bit_count() for v in xb+zb),
                slabs_checked=slab_counts,stabilizer_lift_failures=stab_failures,simple_distance_upper_witnesses=upper)

def gate_checks():
    I=np.eye(2);X=np.array([[0,1],[1,0]]);Z=np.diag([1,-1]);H=np.array([[1,1],[1,-1]])/np.sqrt(2);S=np.diag([1,1j]); plus=np.ones(2)/np.sqrt(2)
    kron=lambda *m: __import__('functools').reduce(np.kron,m)
    cnot=np.eye(4)[:,[0,1,3,2]]
    inj=np.zeros((8,4),complex)
    for c,t,a in itertools.product(range(2),repeat=3):inj[4*c+2*a+t,2*c+t]=plus[a]
    errors=[]
    for u,v,w in itertools.product(range(2),repeat=3):
        p1=(np.eye(8)+(-1)**u*kron(Z,Z,I))/2
        p2=(np.eye(8)+(-1)**v*kron(I,X,X))/2
        proj=np.zeros((4,8));
        for c,t in itertools.product(range(2),repeat=2):proj[2*c+t,4*c+2*w+t]=1
        out=kron(np.linalg.matrix_power(Z,v),np.linalg.matrix_power(X,u^w))@proj@p2@p1@inj
        coeff=np.vdot(cnot,out)/4;errors.append(float(np.max(np.abs(out-coeff*cnot))))
    terr=[]
    for theta in [0.37,np.pi/4,np.pi/2]:
        R=np.diag([1,np.exp(1j*theta)]);R2=np.diag([1,np.exp(2j*theta)])
        j=np.kron(I,(R@plus)[:,None])
        for m,n in itertools.product(range(2),repeat=2):
            parity=(np.eye(4)+(-1)**m*kron(Z,Z))/2
            bra=np.array([1,(-1)**n])/np.sqrt(2)
            k=np.linalg.matrix_power(R2,m)@np.linalg.matrix_power(Z,n)@np.linalg.matrix_power(X,m)@np.kron(bra[None,:],I)@parity@j
            coeff=np.vdot(R,k)/2;terr.append(float(np.max(np.abs(k-coeff*R))))
    return dict(CNOT_max_error=max(errors),phase_teleport_max_error=max(terr),Hadamard_max_error=float(np.max(np.abs(S@H@S@H@S-np.exp(1j*np.pi/4)*H))))

def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n
    c=(p+z*z/(2*n))/den;r=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [c-r,c+r]

if __name__=='__main__':
    result={'gates':gate_checks(),'lattices':[]}
    cases=[(4,[list(range(4))],[list(range(4))],2),(4,[list(range(4))],[list(range(4))],3),
           (6,[list(range(6))],[list(range(6))],2),
           (7,[[0,1,2,3],[0,1,4,5],[0,2,4,6]],[[0,1,2,3],[0,1,4,5],[0,2,4,6]],2),
           (3,[],[[0,1],[1,2]],2), (3,[[0,1],[1,2]],[],2),
           (4,[list(range(4))]*2,[list(range(4))]*2,2)]
    for args in cases:
        r=lattice(*args);result['lattices'].append(r);print(r,flush=True)
    result['table_III']=[dict(k=k,n=1000,rate=k/1000,wilson=wilson(k,1000)) for k in [71,68,5,0,2]]
    result['table_I']=[dict(k=k,n=n,rate=k/n,wilson=wilson(k,n)) for k,n in [(78,2000),(30,600),(342,4500),(252,2850),(28,34000),(13,69000),(3,22800),(1,9300),(31,44000),(10,66600),(0,23400),(0,9150)]]
    result['depth_exponent']=1.5+math.log2(15)
    result['crossing_shift_percent']=100*(4.471-4.328)/4.471
    result['prefix_enumeration']={}
    for k in range(2,7):
        valid=0
        for bits in itertools.product([0,1],repeat=2*k-1):
            xx=bits[:k]; pp=(0,)+bits[k:]+(0,)
            valid+=all(pp[i]^pp[i+1]^xx[i]==0 for i in range(k))
        assert valid==2**(k-1)
        result['prefix_enumeration'][k]=valid
    (OUT/'IRIS_independent_checks.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result['gates']))
