"""Exact Gaussian-rational helical and scaling controls; no PDE simulation."""
from fractions import Fraction as F
from pathlib import Path
import json


def z(a=0,b=0):return (F(a),F(b))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[0],-a[1])
def sub(a,b):return add(a,neg(b))
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,s):return (a[0]*s,a[1]*s)
def dot(a,b):
    r=z()
    for x,y in zip(a,b):r=add(r,mul(x,y))
    return r

def cross(a,b):
    return [sub(mul(a[1],b[2]),mul(a[2],b[1])),sub(mul(a[2],b[0]),mul(a[0],b[2])),sub(mul(a[0],b[1]),mul(a[1],b[0]))]

def curl(k,v):return [mul(z(0,1),x) for x in cross(k,v)]
def realvec(*xs):return [z(x) for x in xs]
def serial(v):return [[str(a),str(b)] for a,b in v]


def run():
    # Unnormalized vectors are sqrt(2) times the parent54 unit inputs.
    k,p=realvec(1,0,0),realvec(0,2,0)
    h,g=[z(),z(1),z(0,1)],[z(0,1),z(),z(1)]
    assert dot(k,h)==z() and dot(p,g)==z()
    assert curl(k,h)==h and curl(p,g)==[scale(x,2) for x in g]
    q=[add(x,y) for x,y in zip(k,p)]
    n=[mul(z(0,1),add(mul(dot(p,h),gg),mul(dot(k,g),hh))) for hh,gg in zip(h,g)]
    assert n==[z(-2),z(-1),z(0,1)]
    qnorm2=dot(q,q);assert qnorm2==z(5)
    v=[sub(nn,scale(mul(qq,dot(q,n)),F(1,5))) for nn,qq in zip(n,q)]
    assert v==[z(F(-6,5)),z(F(3,5)),z(0,1)] and dot(q,v)==z()
    cv=curl(q,v);assert cv==[z(-2),z(1),z(0,3)]
    # Pure positive output would require C_q v=sqrt(5)v; third component
    # would require 3=sqrt(5), contradicted by exact squared integers.
    assert 3**2!=5
    production=3*(F(3,2)+1)-3
    hyper_rhs=1+F(5,4)+F(9,4)
    ordinary_rhs=1+1+2
    assert production==F(9,2)==hyper_rhs and ordinary_rhs==4
    s0six=F(16,3)**3
    matrix_constant=F(27,4)*s0six;assert matrix_constant==1024
    return {"status":"PASS","scope":"Exact datum-level Fourier algebra and scaling; not a PDE simulation, actual trajectory, source payment, or global closure.","normalization":"physical k; input vectors unnormalized, sqrt(2) times parent54 inputs", "convection":serial(n),"leray_output":serial(v),"curl_output":serial(cv),"opposite_sector_nonzero":"pure positive would imply 9=5", "fixed_energy_scaling":{"P":str(production),"beta_5_over_4_rhs":str(hyper_rhs),"beta_1_same_format_rhs":str(ordinary_rhs)},"matrix_G_minus_I_coefficient":str(matrix_constant)}

if __name__=="__main__":
    result=run()
    out=Path(__file__).resolve().parents[2]/"artifacts/unconditional-navier/helical-controls.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
