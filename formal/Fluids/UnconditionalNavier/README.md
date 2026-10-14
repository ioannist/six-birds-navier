# Independent formal controls

Absorption.lean proves the exact real scalar orthogonal absorption law.
It supplies no Navier source premise, no PDE adapter and no unconditional
root. It does not import or alter the older Fluids scaffold.

Run from formal/: `lake env lean Fluids/UnconditionalNavier/Absorption.lean`.
Any future physical source bridge must preserve the original target and
report its actual transitive hypotheses separately from a successful build.

StrainSign.lean proves the scalar ordered trace-free triple inequality
`-4*a*b*c <= 2*max b 0*(a^2+b^2+c^2)`. Run from formal/:
`lake env lean Fluids/UnconditionalNavier/StrainSign.lean`.
The matrix eigenvalue substitution, actual-path Piola/Betchov identity and
conditional dose continuation are written research proofs, outside this
declaration. No terminal strain-dose supplier or full PDE theorem is checked
by this scalar file.
