#!/usr/bin/env python3
"""Exact exponent/geometry controls; no original trajectory is simulated."""
from fractions import Fraction as F
import json
from pathlib import Path

def controls():
 family={}
 for q in (F(2),F(5,2),F(3),F(7,2),F(4)):
  aq=F(3,2)-3/q;bq=-F(1,4)+3/q
  assert aq+bq==F(5,4) and bq/(1-aq)>1
  family[str(q)]={'a':str(aq),'b':str(bq),'holder_kernel_power':str(bq/(1-aq)),
                  'fixed_native_budget_band_majorant_power':str(1-aq-bq)}
 # Low frequency a=(1,0,0); b=(R,1,0); c=-a-b.
 geometry={}
 for R in (4,16,64):
  a=(1,0,0);b=(R,1,0);c=tuple(-x-y for x,y in zip(a,b))
  short=sum((x+y)**2 for x,y in zip(b,c))
  long=sum((x-y)**2 for x,y in zip(b,c))
  assert short==1 and long==(2*R+1)**2+4
  geometry[str(R)]={'reflected_segment_length_squared':short,'unreflected_segment_length_squared':long}
 # For gamma3: p=k+1, normalized objective 2d+(2/k)d^-k.
 # Its derivative is 2-2/d^p, so the root has d^p=1.
 k=F(13,4);p=F(17,4);wrong_z=k/p
 wrong_derivative=2-2/wrong_z
 assert wrong_derivative==F(-8,13)
 return {'scope':'Exact rational exponent, triad-geometry and minimizer controls; no PDE event realization or signed source bound.',
         'native_interpolation_family':family,'triad_geometry':geometry,
         'gamma3_split_minimizer':{'correct_normalized_d_power':str(F(1)),
             'incorrect_old_d_power':str(wrong_z),'objective_derivative_when_d_power_equals_old_value':str(wrong_derivative)},
         'unpaid':'Original signed cone-tip source, actual cadence/disjointness and pressure/critical supplier.'}
if __name__=='__main__':
 result=controls()
 target=Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/source-window-controls.json'
 target.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
