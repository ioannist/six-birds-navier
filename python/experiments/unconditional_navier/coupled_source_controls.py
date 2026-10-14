"""Exact coupled-source and tail controls; no PDE simulation or source payment."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as s


def controls():
    a,b,c,d,e,f=s.symbols("a b c d e f")
    matrix=s.Matrix([[a,b,c],[b,d,e],[c,e,f]])
    cofactor=2*matrix**2-2*s.trace(matrix)*matrix+(s.trace(matrix)**2-s.trace(matrix**2))*s.eye(3)
    assert all(s.expand(v)==0 for v in cofactor-2*matrix.adjugate())
    positive=s.diag(0,1,1).adjugate()[0,0]
    negative=s.diag(0,1,-1).adjugate()[0,0]
    assert positive==1 and negative==-1
    entries=s.symbols("a0:8");A=s.Matrix([[entries[0],entries[1],entries[2]],[entries[3],entries[4],entries[5]],[entries[6],entries[7],-entries[0]-entries[4]]])
    T=s.Matrix(3,3,s.symbols("t0:9"))
    def axial(M):
        return s.Matrix([sum(s.LeviCivita(i,j,k)*M[k,j] for j in range(3) for k in range(3)) for i in range(3)])
    assert all(s.expand(v)==0 for v in axial(-A*T-T*A.T)-A.T*axial(T))
    x,y,z=s.symbols("x y z");coords=(x,y,z)
    def curl(v):
        return s.Matrix([s.diff(v[2],y)-s.diff(v[1],z),s.diff(v[0],z)-s.diff(v[2],x),s.diff(v[1],x)-s.diff(v[0],y)])
    u=curl(s.Matrix([0,0,-x*y*y/2]));omega=curl(u)
    assert sum(s.diff(u[k],coords[k]) for k in range(3))==0
    rho_t=-s.Matrix([sum(u[k]*s.diff(omega[i],coords[k]) for k in range(3)) for i in range(3)])
    assert rho_t==s.Matrix([0,0,x*y]) and curl(rho_t)==s.Matrix([x,-y,0])
    assert F(27,32)*F(16,3)**3==128
    exponents=[]
    for p in (F(9,2),F(6)):
        near=F(3,2)-6/p
        far=-6/p
        flux_near=-1+near+3*(F(1,2)-1/p)
        flux_far=-1+far+3*(1-1/p)
        assert flux_near==flux_far==2-9/p
        exponents.append(dict(p=str(p),pressure_near_L2_radius_power=str(near),pressure_far_radius_power=str(far),full_flux_radius_power=str(flux_near),diffusion_cutoff_radius_power=str(1-6/p)))
    # Disjoint annular norm profiles W_r=r^(-7/12)F(x/r), r=2^(12n).
    # They are not stationary solutions or actual record limits.
    profiles=[dict(n=n,radius=2**(12*n),Z_factor=str(F(1,2**(2*n))),L6_sixth_factor=str(F(1,2**(6*n))),absolute_cubic_cutoff_price_factor=2**(3*n)) for n in (1,2,3,4)]
    tau,R,nu=s.symbols("tau R nu",positive=True)
    absolute=s.integrate((nu*tau)**(-s.Rational(5,4)),(tau,R,s.oo))
    signed_dissipation=nu*(4*s.pi*nu)**(-s.Rational(3,2))*s.integrate(tau**(-s.Rational(3,2)),(tau,R,s.oo))
    assert s.simplify(absolute-4/nu*(nu*R)**(-s.Rational(1,4)))==0
    assert s.simplify(signed_dissipation-1/(4*s.pi**s.Rational(3,2)*s.sqrt(nu*R)))==0
    return dict(status="PASS",scope="Exact polynomial identities, initial inner-field source creation, scaling and norm-profile controls only; no PDE simulation, stationary realization, endpoint source or all-data closure.",cofactor_sign_controls=[int(positive),int(negative)],axial_reaction_residual="0",actual_inner_polynomial_rho_t=[str(v) for v in rho_t],local_Young_coefficient=128,stationary_flux_exponents=exponents,annular_norm_profiles=profiles,absolute_heat_tail=str(absolute),signed_heat_dissipation_tail=str(signed_dissipation))


if __name__=="__main__":
    result=controls()
    target=Path(__file__).resolve().parents[2]/"artifacts/unconditional-navier/coupled-source-controls.json"
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(status=result["status"],scope=result["scope"],artifact=str(target)),indent=2))
