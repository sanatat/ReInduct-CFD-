# Governing equations, coupling decisions and limitations

## Existing C++ reduced-order model
- Gas species: epsilon*dc_i/dt + d(U*c_i)/dz = nu_i*r.
- Solid conversion: nNiO*dX/dt = r; 0 <= X <= 1.
- Gas energy: epsilon*Cp_vol_g*dTg/dt + U*Cp_vol_g*dTg/dz = hvol*(Ts-Tg).
- Solid energy: Cp_vol_s*dTs/dt = k*d²Ts/dz² + q_ind(z) - hvol*(Ts-Tg) - hLoss*(Ts-Tambient) - ΔH*r.
- Illustrative rate r = nNiO*A*exp(-Ea/(R*Ts))*(cH2/cRef)*(1-X).
- NOTE: the prescribed `q_ind(z)` in C++ is not linked to magnetic solution; do not claim otherwise.

## New independently verified induction benchmark
- Harmonic magnetic diffusion for 1-D conductive slab: H'' = i*omega*mu*sigma*H.
- BCs H(0)=H0, H'(L)=0; exact complex field H(x) = H0*cosh(k*(L-x))/cosh(k*L), k=sqrt(i*omega*mu*sigma)=(1+i)/delta.
- delta = sqrt(2/(omega*mu*sigma)).
- Complex current amplitude J = -dH/dx; cycle-averaged Joule heat density q=|J|²/(2*sigma).
- Model is magnetoquasistatic and neglects displacement current, temperature dependence of sigma, coil geometry, and Lorentz feedback.

## New sequential one-way thermal model
- rhoCp*dT/dt = kthermal*d²T/dx² + q_EM(x) - hvol*(T-Tgas).
- No-flux x boundaries. Implicit finite-volume update (tridiagonal solve).
- Prescribed constant Tgas, not a coupled CFD flow field. Coupling only in EM -> thermal direction.

## Full target equations (NOT YET IMPLEMENTED)
- Continuity: d(rho*epsilon)/dt + div(rho*U) = 0.
- Momentum: d(rho*U)/dt + div(rho*U*U) = -grad(p) + div(tau) - Darcy/Forchheimer resistance + JxB.
- Species: d(epsilon*rho*Yk)/dt + div(rho*U*Yk) = div(epsilon*rho*Dk*grad(Yk)) + omega_dot_k.
- Energy: d(rho*h)/dt + div(rho*U*h) = div(k*grad(T)) + qJoule + qreaction + interfacial transfer + radiation.
- Electromagnetism: curl(E)=-dB/dt; curl(H)=J; B=mu*H; J=sigma*(E+U x B).
- Surface solid reactions: species stoichiometric conservation plus reversible kinetics (fitted to experimental evidence).

## Validation ladder
1. Unit checks, positivity, exact manufactured or analytical solutions.
2. Closed inventory and flux budgets.
3. Mesh, timestep, nonlinear residual independence.
4. Benchmarks (cavity, heated slab, reacting plug flow).
5. Published experimental validation for NiO reduction kinetics.
6. Battery black-mass dataset validation before industrial claims.
