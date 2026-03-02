# Fluids conventions (NS)

## Domain and Fourier conventions
- Domain: T^d with d in {2,3}
- Fourier series: u(x) = sum_k u_hat(k) * exp(i k·x), k in Z^d
- Discrete grid: use numpy.fft.fftfreq with k = 2*pi*fftfreq(N, d=L/N)

## Divergence-free constraint and Leray projection
- For k != 0: P_k = I - (k k^T)/|k|^2
- For k = 0: P_0 = I (leave the zero mode unchanged)

## Operators
- Laplacian multiplier: Delta -> -|k|^2
- Stokes operator: A = -P Delta, so A u_hat = |k|^2 P_k u_hat
- Hyperviscous multiplier: (-Delta)^alpha -> |k|^(2 alpha)

## Bilinear term notation
- B(u,v) = P(u · grad v)
