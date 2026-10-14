#!/usr/bin/env python3
"""Exact periodic seed controls and scalar constants; no nonlinear PDE simulation."""
from fractions import Fraction as F
from pathlib import Path
import json

ZERO=(F(0),F(0)); ONE=(F(1),F(0))
def ca(a,b):return a[0]+b[0],a[1]+b[1]
def cm(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def add(a,b):
 c=dict(a)
 for k,v in b.items():c[k]=ca(c.get(k,ZERO),v)
 return {k:v for k,v in c.items() if v!=ZERO}
def scale(a,q):return {k:cm(v,(F(q),F(0))) for k,v in a.items() if v!=ZERO}
def mul(a,b):
 c={}
 for k,v in a.items():
  for l,w in b.items():
   z=tuple(x+y for x,y in zip(k,l));c[z]=ca(c.get(z,ZERO),cm(v,w))
 return {k:v for k,v in c.items() if v!=ZERO}
def deriv(a,j):return {k:cm(v,(F(0),F(k[j]))) for k,v in a.items() if k[j]}
def cos(k):return {tuple(k):(F(1,2),F(0)),tuple(-x for x in k):(F(1,2),F(0))}
def sin(k):return {tuple(k):(F(0),F(-1,2)),tuple(-x for x in k):(F(0),F(1,2))}
def curl(u):return [add(deriv(u[2],1),scale(deriv(u[1],2),-1)),add(deriv(u[0],2),scale(deriv(u[2],0),-1)),add(deriv(u[1],0),scale(deriv(u[0],1),-1))]
def production(u):
 w=curl(u);s=[[scale(add(deriv(u[i],j),deriv(u[j],i)),F(1,2)) for j in range(3)] for i in range(3)]
 p={}
 for i in range(3):
  for j in range(3):p=add(p,mul(mul(w[i],s[i][j]),w[j]))
 return p.get((0,0,0),ZERO)

def main():
 x=(1,0,0);y=(0,1,0);xy=(1,1,0)
 u=[sin(y),{},add(sin(x),scale(sin(xy),-1))]
 a=[{},add(scale(cos(x),-1),cos(xy)),scale(cos(y),-1)]
 assert curl(a)==u
 assert not add(add(deriv(u[0],0),deriv(u[1],1)),deriv(u[2],2))
 mean=production(u);assert mean==(F(1,4),F(0))
 assert production([scale(v,-1) for v in u])==(F(-1,4),F(0))
 c04=F(4)**6/F(3)**3/4;cg=F(27,128)*c04;assert cg==8
 high=F(27,2048)*4**4*c04;assert high==128
 result={'scope':'Exact periodic algebra plus whole-space compact-curl averaging proof in parent47. Not a PDE solution, endpoint estimate or closure.', 'periodic_curl_and_divergence':'exactly verified','periodic_mean_production':str(mean[0]),'phase_reversal_mean':'-1/4','C0_fourth_power':str(c04),'native_growth_coefficient':str(cg),'split_high_frequency_coefficient':str(high),'event_A_fourth_power':str(4**4*c04*2**5),'profile_lower_coefficient_over_pi':'19/30','unresolved':['actual signed weighted source','all-event or long-duration continuation return']}
 out=Path('python/artifacts/unconditional-navier/weighted-velocity-controls.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
