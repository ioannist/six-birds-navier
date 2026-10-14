#!/usr/bin/env python3
"""Exact actual curl-Gaussian integrated bending controls; requires SymPy.

Finite initial-state controls only; not a universal integrated sign theorem.
The Gaussian moment algorithm is justified in round09/parent-controls.md.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import sympy as s

x=s.symbols('x y z',real=True)
def derivative(p,j,rate=1):
    return s.expand(s.diff(p,x[j])-2*rate*x[j]*p)
def curl(p):
    return s.Matrix([derivative(p[2],1)-derivative(p[1],2),
                     derivative(p[0],2)-derivative(p[2],0),
                     derivative(p[1],0)-derivative(p[0],1)])
def normalized_moment(p):
    # All tested integrated quartic fields have rate4.
    total=s.Integer(0)
    for powers,c in s.Poly(s.expand(p),*x).terms():
        if any(n%2 for n in powers):
            continue
        total+=c*s.prod(s.factorial2(n-1)/8**(n//2) for n in powers)/8
    return s.factor(total)
def coefficients(u):
    assert s.expand(sum(derivative(u[j],j) for j in range(3)))==0
    w=curl(u)
    v=[s.Matrix([derivative(p,j) for p in w]) for j in range(3)]
    J=[w.cross(a).applyfunc(s.expand) for a in v]
    a=s.Integer(0);b=s.Integer(0)
    for j in range(3):
        dw=v[j]
        for k in range(3):
            dJ=s.Matrix([derivative(p,j,2) for p in J[k]])
            ddw=s.Matrix([derivative(p,j) for p in v[k]])
            a+=normalized_moment(dJ.dot(dJ))
            b+=normalized_moment(J[k].dot(dw.cross(ddw)))
    return tuple(s.factor(z) for z in (a,b,-2*a-4*b))

U=s.Matrix([4*x[0]*x[2]**2-2*x[0],-2*x[1]*x[2]**2+x[1],
            -4*x[0]**2*x[2]+2*x[1]**2*x[2]+x[2]])
states=[('original-U',U,(s.Rational(8628471,16384),s.Rational(10183671,65536),s.Rational(-27440613,16384))),
        ('potential-(0,0,1)',curl(s.Matrix([0,0,1])),(s.Rational(627,4),s.Rational(627,16),s.Rational(-1881,4))),
        ('potential-(x²,0,0)',curl(s.Matrix([x[0]**2,0,0])),(s.Rational(293395,16384),s.Rational(216595,65536),s.Rational(-803385,16384))),
        ('potential-(yz,x,1)',curl(s.Matrix([x[1]*x[2],x[0],1])),(s.Rational(34353515,16384),s.Rational(59050971,65536),s.Rational(-127758001,16384)))]
result={}
for name,u,expected in states:
    found=coefficients(u)
    assert found==expected,(name,found,expected)
    assert found[0]>=0 and found[2]<0
    result[name]={key:str(value) for key,value in zip(('A/pi^(3/2)','B/pi^(3/2)','V/pi^(3/2)'),found)}
    print(name,result[name],flush=True)

# Independently evaluate the normalized density derivative at the original
# compactifiable solenoidal shear jet. These are local actual initial controls.
w=s.Matrix([1+x[0]**2/4,x[0]-x[0]*x[1]/2,0])
v=[w.diff(z) for z in x]
J=[w.cross(a) for a in v]
K=s.expand(sum(a.dot(a) for a in J));L2=s.expand(w.dot(w))
at0={z:0 for z in x}
lapw=s.Matrix([sum(s.diff(p,z,2) for z in x) for p in w])
dwJ=[lapw.cross(v[j])+w.cross(lapw.diff(x[j])) for j in range(3)]
heatK=s.expand(2*sum(J[j].dot(dwJ[j]) for j in range(3)))
assert K.subs(at0)==1 and sum(s.diff(K,z,2) for z in x).subs(at0)==0
assert heatK.subs(at0)==1 and w.dot(lapw).subs(at0)==s.Rational(1,2)
normalizations={}
for p in (1,2,4):
    G=L2**s.Rational(p-4,2)*K
    heatG=s.Rational(p-4)*L2**s.Rational(p-6,2)*w.dot(lapw)*K+L2**s.Rational(p-4,2)*heatK
    lapG=sum(s.diff(G,z,2) for z in x)
    source=s.simplify((heatG-lapG).subs(at0))
    assert source==5-p
    normalizations[str(p)]={'G(0)':str(G.subs(at0)),'D_vis G /nu':str(heatG.subs(at0)),
                            'DeltaG(0)':str(lapG.subs(at0)),
                            '(D-nuDelta)G_vis/nu':str(source)}
artifact={'integrated_initial_controls':result,'normalized_shear_controls':normalizations,
          'sympy_version':s.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'scope':'Finite actual solenoidal Schwartz Gaussian controls and compactifiable shear jet only; no universal integrated sign, trajectory sign or root closure.'}
out=Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/bending-integrated-controls.json'
out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(artifact,indent=2)+'\n')
print(json.dumps(normalizations,indent=2))
