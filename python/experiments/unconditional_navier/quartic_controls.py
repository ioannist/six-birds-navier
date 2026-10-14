#!/usr/bin/env python3
"""Fresh ONE-heated-receiver original quartic; exact full Leray projection.

No execution/import of the historical all-heated-memory script.
Angular-frequency transform exp(-ik.x), Parseval (2pi)^-3.
"""
from functools import lru_cache
from pathlib import Path
from math import isqrt
import sympy as s
import json,hashlib
OUT=Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/quartic-controls.json'
x=s.symbols('x y z',real=True); k=s.symbols('k1 k2 k3',real=True);r=s.symbols('r',positive=True)
U=s.Matrix([4*x[0]*x[2]**2-2*x[0],-2*x[1]*x[2]**2+x[1],-4*x[0]**2*x[2]+2*x[1]**2*x[2]+x[2]])
def d(p,j,a=1):return s.expand(s.diff(p,x[j])-2*a*x[j]*p)
B=s.Matrix([s.expand(sum(U[j]*d(U[i],j) for j in range(3))) for i in range(3)])
assert s.expand(sum(d(U[j],j) for j in range(3)))==0
@lru_cache(None)
def mon(n,a):
 if sum(n)==0:return s.Integer(1)
 j=next(j for j,v in enumerate(n) if v);m=list(n);m[j]-=1;q=mon(tuple(m),a)
 return s.expand(s.I*(s.diff(q,k[j])-k[j]*q/(2*a)))
def ft(p,a):return s.expand(sum(c*mon(n,s.Rational(a)) for n,c in s.Poly(p,*x).terms()))
def angular(p):
 total=0
 for n,c in s.Poly(s.expand(p),*k).terms():
  if any(v%2 for v in n):continue
  total+=c*s.prod(s.factorial2(v-1) for v in n)/s.factorial2(sum(n)+1)*r**sum(n)
 return s.factor(total)
def radial(mean,alpha):
 expr=s.cancel(r**2*mean)
 assert not s.denom(expr).has(r)
 value=0
 for (n,),c in s.Poly(s.expand(expr),r).terms():
  assert n>=0
  h=s.Rational(n+1,2);value+=c*s.gamma(h)/(2*alpha**h)
 return s.simplify(value)
def average(p,a):
 return s.factor(sum(c*s.prod(s.factorial2(v-1)/(2*s.Rational(a))**(v//2) for v in n)
    for n,c in s.Poly(s.expand(p),*x).terms() if not any(v%2 for v in n)))
def heat_poly(p,a,t):
 denom=1+4*a*t;center=1/denom;variance=2*t/denom;out=0
 for n,c in s.Poly(p,*x).terms():
  factors=[]
  for j,m in enumerate(n):
   factors.append(sum(s.binomial(m,2*ell)*s.factorial2(2*ell-1)*variance**ell*(center*x[j])**(m-2*ell)
                      for ell in range(m//2+1)))
  out+=c*s.prod(factors)
 return s.expand(out),denom**s.Rational(-3,2),a/denom
bh=s.Matrix([ft(p,2) for p in B])
def pairing(left,a,heat_extra=0):
 """Pair (Gaussian left polynomial,ratea) against FULL P B, plus heat factor."""
 lh=s.Matrix([ft(p,a) for p in left])
 raw=angular(sum(s.conjugate(bh[i])*lh[i] for i in range(3)))
 long=angular(sum(k[i]*s.conjugate(bh[i]) for i in range(3))*sum(k[i]*lh[i] for i in range(3)))/r**2
 pref=s.simplify((s.pi/a)**s.Rational(3,2)*(s.pi/2)**s.Rational(3,2)*4*s.pi/(2*s.pi)**3)
 alpha=1/(4*a)+s.Rational(1,8)+heat_extra
 rawvalue=s.simplify(pref*radial(raw,alpha))
 longvalue=s.simplify(pref*radial(long,alpha))
 return rawvalue,longvalue,s.simplify(rawvalue-longvalue),{'prefactor':str(pref),'radial_alpha':str(alpha),'mean_raw':str(raw),'mean_longitudinal':str(s.factor(long))}

def compute(beta):
 lam=s.sqrt(1/(1+8*beta));a=lam**2
 hu=s.Matrix([s.expand(lam**6*p.subs({x[j]:lam*x[j] for j in range(3)},simultaneous=True)) for p in U])
 mixed=s.Matrix([s.expand(sum(d(U[h],i)*hu[h] for h in range(3))-sum(U[j]*d(hu[i],j,a) for j in range(3))) for i in range(3)])
 first=pairing(B,s.Integer(2),2*beta)
 second=pairing(mixed,1+a)
 # Independent physical unprojected check of both factors.
 hbp=[]
 for p in B:
  q,amp,rate=heat_poly(p,s.Integer(2),2*beta);hbp.append(q)
 physical_first=s.simplify(amp*average(sum(hbp[i]*B[i] for i in range(3)),rate+2)*(s.pi/(rate+2))**s.Rational(3,2))
 physical_second=s.simplify(average(mixed.dot(B),3+a)*(s.pi/(3+a))**s.Rational(3,2))
 assert s.simplify(first[0]-physical_first)==0
 assert s.simplify(second[0]-physical_second)==0
 q=s.simplify(-first[2]-second[2])
 return {'beta':str(beta),'lambda':str(lam),'receiver_rate':str(a),'mixed_rate':str(1+a),
   'heated_B':{'raw':str(first[0]),'longitudinal':str(first[1]),'projected':str(first[2]),**first[3]},
   'mixed_E':{'raw':str(second[0]),'longitudinal':str(second[1]),'projected':str(second[2]),**second[3]},
   'Q_beta':str(q),'Q_beta_over_pi_3half':str(s.simplify(q/s.pi**s.Rational(3,2))),
   'physical_unprojected_checks':'Both physical polynomial heat/Gaussian pairings match exact Fourier raw pairings.'},q
print('Computingbeta3/8withfullprojection',flush=True)
res,q=compute(s.Rational(3,8));print('Q_beta='+str(q),flush=True)
print('Checking beta0 exactcancellation',flush=True)
zero,qzero=compute(s.Integer(0));assert qzero==0
coefficient=s.expand(s.simplify(q/s.pi**s.Rational(3,2)))
# Exact radical enclosure; never use a decimal to certify sign.
rad=s.sqrt(13);c=coefficient.coeff(rad);a=s.simplify(coefficient-c*rad)
assert a.is_Rational and c.is_Rational
scale=10**30;n=isqrt(13*scale*scale);lo=s.Rational(n,scale);hi=s.Rational(n+1,scale)
assert lo**2<13<hi**2
ql=s.simplify(a+c*(lo if c>=0 else hi));qh=s.simplify(a+c*(hi if c>=0 else lo))
sign='positive' if ql>0 else 'negative' if qh<0 else 'undetermined'
res.update(beta0_control=zero,radical_enclosure=[str(ql),str(qh)],sign=sign,
  amplitude_parity='G and F each quadratic, Q(bU)=b^4Q(U), identical for +/-U.',
  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=s.__version__,
  scope='Exact original initial-state ONEheatedreceiver quartic only; not trajectorydrift, infiniteheatmemory, all-data source or mechanized PDE theorem.')
assert sign=='negative'
assert c>0 and a<0 and s.simplify(a*a-13*c*c)>0
res['squared_radical_sign_certificate']=str(s.factor(a*a-13*c*c))
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
