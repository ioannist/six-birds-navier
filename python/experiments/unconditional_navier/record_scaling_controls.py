"""Exact norm controls for record scaling and annular budgets; no PDE simulation."""
from fractions import Fraction as F
from pathlib import Path
import json

def controls():
    rows=[]
    for s in (2,4,8,16,32):
        # L=s^2; u0=L^(3/2)W(Lx), Z_W=1. Coefficients multiply base norms.
        b=s**4
        q=F(1,s**4)
        p=F(1,s**3)
        critical=s
        # alpha=3/4, tau=s^-4, B=s^3: normalized backward Z ratio fourth power.
        slow_ratio_fourth=(F(1)+F(1,s**2))**-3
        rows.append(dict(s=s,physical_energy_factor=1,physical_Z=b,
                         normalized_Z=1,normalized_Q_factor=str(q),
                         normalized_P_factor=str(p),normalized_L3_norm_factor=critical,
                         slow_scalar_record_ratio_fourth=str(slow_ratio_fourth)))
    annuli=[]
    for n in (1,2,4,8,16):
        z=sum((F(1,4**j) for j in range(1,n+1)),F(0))
        d=sum((F(1,4**(3*j)) for j in range(1,n+1)),F(0))
        assert z < F(1,3) and d < F(1,63)
        annuli.append(dict(count=n,critical_L3_mass_factor=n,enstrophy_factor=str(z),
                           laplacian_and_L6_mass_factor=str(d)))
    return dict(status="PASS",scope="Exact symbolic scaling and disjoint profile controls; no actual singular path, same-datum no-escape conclusion, source payment, or NS simulation.",
                record_data_family=rows,annular_profiles=annuli,
                kernel_scaling=dict(K_squared_L2="lambda^3",normalized_raw="lambda^2",
                                    J_squared_L2="lambda",normalized_covariant="1"))

if __name__ == "__main__":
    result=controls()
    target=Path(__file__).resolve().parents[2]/"artifacts/unconditional-navier/record-scaling-controls.json"
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(status=result["status"],scope=result["scope"],artifact=str(target)),indent=2))
