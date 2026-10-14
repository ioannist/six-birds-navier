#!/usr/bin/env python3
"""Exact original initial-state controls; no trajectory/source estimate."""
from fractions import Fraction as F
from pathlib import Path
from itertools import permutations
import json
from heat_ledger_controls import add,scale,multiply,derivative,gaussian_average

def lap(poly,rate=1):
 return add(*(derivative(derivative(poly,j,rate),j,rate) for j in range(3)))
def curl(vector,rate=1):
 return [add(derivative(vector[(i+2)%3],(i+1)%3,rate),
             scale(derivative(vector[(i+1)%3],(i+2)%3,rate),-1)) for i in range(3)]
def advect(a,b,rate=1):
 return [add(*(multiply(a[j],derivative(b[i],j,rate)) for j in range(3))) for i in range(3)]
def pair(a,b,rate=3):
 return gaussian_average(add(*(multiply(x,y) for x,y in zip(a,b))),rate)
def controls():
 u=[{(1,0,2):F(4),(1,0,0):F(-2)},
    {(0,1,2):F(-2),(0,1,0):F(1)},
    {(2,0,1):F(-4),(0,2,1):F(2),(0,0,1):F(1)}]
 assert not add(*(derivative(u[i],i) for i in range(3)))
 w=curl(u);lw=[lap(v) for v in w]
 transport=-pair(lw,[lap(v,2) for v in advect(u,w)])
 stretching=pair(lw,[lap(v,2) for v in advect(w,u)])
 B=advect(u,u);curlB=curl(B,2)
 total=-pair(lw,[lap(v,2) for v in curlB])
 lap3=[lap(lap(lap(v))) for v in u]
 assert transport==F(23936,81) and stretching==F(35200,81)
 assert total==transport+stretching==F(19712,27)==pair(lap3,B)
 # Actual Schwartz pure swirl datum. Its quadratic convection is formed
 # from original u; derivative rate2 is essential for curl(B).
 swirl=[{(0,1,0):F(-2)},{(1,0,0):F(2)},{}]
 assert not add(*(derivative(swirl[i],i) for i in range(3)))
 swirl_B=advect(swirl,swirl)
 assert swirl_B==[{(1,0,0):F(-4)},{(0,1,0):F(-4)},{}]
 swirl_curlB=curl(swirl_B,2)
 assert swirl_curlB==[{(0,1,1):F(-16)},{(1,0,1):F(16)},{}]
 # S=[[-1,1,0],[1,0,0],[0,0,1]] has characteristic polynomial
 # (lambda-1)(lambda^2+lambda-1)=lambda^3-2lambda+1.
 strain=[[-1,1,0],[1,0,0],[0,0,1]]
 def pmul(a,b):
  out=[F(0)]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):out[i+j]+=x*y
  return out
 chi=[F(0)]*4
 for perm in permutations(range(3)):
  sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
  term=[F(sign)]
  for i,j in enumerate(perm):term=pmul(term,[F(-strain[i][j]),F(int(i==j))])
  for j,v in enumerate(term):chi[j]+=v
 assert chi==[F(1),F(-2),F(0),F(1)] and sum(chi)==0
 return {'scope':'Exact original Schwartz initial-state nonlinear branches and explicit false-applicability controls, not PDE evolution/global all-data source.',
   'gaussian_coefficients_relative_to_(pi/3)^(3/2)':{
     'transport':str(transport),'stretching':str(stretching),'sum':str(total)},
   'gamma3_weighted_fourier_bridge':'Signed velocity Fourier pairing W=-(T+S)/(2*pi)^6; vorticity nonlinear production is +(T+S)/(2*pi)^6. Viscous dissipation is separate and quadratic.',
   'swirl_control':{'original_convection_polynomial':[['-4x'],['-4y'],['0']],
     'convection_gaussian_rate':2,
     'eta_initial':'0','eta_initial_time_derivative':'-16*z*exp(-2*(r^2+z^2))',
     'point_r1_z1':'-16*exp(-4), nonzero'},
   'compact_stream_swirl_jet':{'strain':[[-1,1,0],[1,0,0],[0,0,1]],
      'c_u_r_over_r':0,'characteristic_polynomial_ascending':[str(v) for v in chi],
      'middle_eigenvalue':'(sqrt(5)-1)/2, strictly between 0 and 1',
      'original_state_bridge':'Parent written cutoff STREAM construction preserves divergence, jet and Schwartz compact support.'},
   'unpaid':'All-data signed cone source and all-data applicability to no-swirl scalar/strain certificate.'}
if __name__=='__main__':
 result=controls()
 artifact=Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/vorticity-controls.json'
 artifact.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
