# Physics equations, implementation and limits — v3

## 1. Navier–Stokes (MAC staggered-grid benchmark)

Constant-density incompressible Navier–Stokes on a 2-D rectangle, x-periodic, no-slip in y:

\[\n\nabla\cdot\mathbf{u}=0, \qquad
\partial_t\mathbf{u}+\mathbf{u}\cdot\nabla\mathbf{u}
=-\nabla(p/\rho)+\nu\nabla^2\mathbf{u}+\mathbf{g}-\boldsymbol{\beta}\mathbf{u}.\n\]

Predict u* with explicit upwind convection and centered viscous terms; solve div(grad p_c)=div(u*)/dt with Fourier/DCT for periodic–Neumann boundary operators; correct u = u* − dt grad(p_c). 2-D x-localized Brinkman resistance is a **prescribed macroscopic resistance**, not resolved solid particles. The `v` component responds to pressure gradients from nonuniform x-resistance.

For the idealized imposed transverse B=(0,B,0), E=0, quasistatic magnetic force is \(f_x/\rho=-\sigma B^2 u/\rho\). The Hartmann profile for this specific idealized assumption is

\[u(y)=\frac{g}{\beta}\left[1-\frac{\cosh(\sqrt{\beta/\nu}(y-H/2))}{\cosh(\sqrt{\beta/\nu}H/2)}\right],\quad\beta=\frac{\sigma B^2}{\rho}.\]

The Hartmann test is independent of the H₂ reactor, whose gas is not assumed to be electrically conducting. **This is not a self-consistent plasma MHD solver.**

## 2. Reactive 2-D finite-volume gas/solid transport

The gas and stationary solid **interpenetrate within the same 2-D cell** with void fraction ε:

\[\epsilon\partial_t C_{H_2}+\nabla\cdot(\mathbf U C_{H_2}-\epsilon D\nabla C_{H_2})=-r,\]
\[\epsilon\partial_t C_{H_2O}+\nabla\cdot(\mathbf U C_{H_2O}-\epsilon D\nabla C_{H_2O})=+r,\]
\[n_{\mathrm{NiO},0}\partial_t X=r.\]

The superficial velocity U comes either from a prescribed plug/Poiseuille field or from the separate 2-D MAC flow computation, passed one-way into a reactor with open inflow/outflow for concentrations. A periodic representative-cell hydrodynamic solution is used as a frozen velocity field; there is no feedback from reactant conversion or heating to momentum.

With effective volumetric gas heat capacity Cp,g and solid heat capacity Cs,bulk:

\[\epsilon C_{p,g}\partial_t T_g + \nabla\cdot(\mathbf U C_{p,g}T_g-k_g\nabla T_g)=h_v(T_s-T_g),\]
\[C_{s,bulk}\partial_t T_s-\nabla\cdot(k_s\nabla T_s)=q_{EM}-\Delta H\,r-h_v(T_s-T_g).\]

Radiation is a boundary heat loss on the upper and lower walls, \(q''=\epsilon_{rad}\sigma_{SB}(T_s^4-T_\infty^4)\). All FV fluxes are constructed once per face, creating an exactly conservative discrete integrated H2/energy bookkeeping in the implemented model. Stability is checked with advective/diffusive CFL limits. The model uses explicit time stepping.

The illustrative rate expression is

\[r=n_0 A e^{-E_a/(RT_s)}(1-X)\max\bigl(0,\,[C_{H_2}-C_{H_2O}/K(T_s)]/C_\mathrm{ref}\bigr),\]

with \(\Delta G(T)=\Delta H(T)-T\Delta S(T)\), \(K(T)=\exp[-\Delta G/(RT)]\), and user-adjustable illustrative thermodynamic parameters. This rate is **not a validated heterogeneous NiO kinetic mechanism**, and the imposed equilibrium-limited cutoff is a simplified driving factor, not thermodynamic phase equilibrium.

## 3. Electromagnetics and temperature feedback

Complex time-harmonic 1-D tangential magnetic field in conducting planar slab:

\[\frac{d}{dy}\left[\frac{1}{\sigma(T)}\frac{dH}{dy}\right]=i\omega\mu H, \quad H(0)=H_0,\quad \partial_y H(L)=0.\]

Calculate \(J=-\partial_y H\), Joule heating \(q_{EM}=|J|^2/(2\sigma)\), and update \(\sigma(T)=\sigma_0/[1+\alpha(T-T_0)]\) using **axially averaged solid temperature** at a prescribed update interval. The heat is projected onto the 2-D grid. This is a reduced sequential EM/T coupling (conductivity changes with T) and ignores actual coil geometry, radiative EM effects, magnetic saturation, mechanical Lorentz feedback and full plasma transport.

## 4. Prescribed-rate radical bookkeeping (separate module)

H₂ → 2H, 2H → H₂ and NiO+2H → Ni+H₂O. Integrates coupled 0-D ODEs with BDF and verifies 
\(2n_{H_2}+n_H+2n_{H_2O}\) and \(n_{H_2O}=n_0X\). All effective rates are illustrative; no electron-impact cross sections, ionization, electron temperature or EEDF.

## 5. Missing realism / next research gates

1. **Real battery black mass:** NMC/LFP/carbon/binder/impurity mixture; multicomponent phases and independent thermodynamic data.
2. **OpenFOAM solver:** compile and run actual reacting species and gas/solid energy models; provide log, `checkMesh`, mesh/refinement tests. The included OpenFOAM cases are only unexecuted scaffolds.
3. **Multiphase:** actual particle trajectories or Eulerian volume fraction momentum equations; separate mesh/interface for multiregion.
4. **3-D coil EM:** realistic induction geometry using a suitable FEM/EM library with volumetric power and Lorentz mapping.
5. **Actual plasma:** electron energy/radiation, chemically validated radical mechanisms and conductivity/transport closure.
6. **Experimental validation:** compare predicted conversion and metal recovery with measured reactor data and identify parameters.

These limits must accompany any public repository, CV, manuscript or research statement.

## 6. Separate-region conjugate heat-transfer benchmark

`multiregion/conjugate1d.py` meshes a low-conductivity layer and a high-conductivity layer with different cell widths and independent interface-adjacent finite-volume cells. It enforces a **single common heat flux** using contact resistance and half-cell conductances:

\[q''=(T_{g,last}-T_{s,first})/[\Delta x_g/(2k_g)+R_c+\Delta x_s/(2k_s)].\]

The steady result is compared against \(q''=(T_L-T_R)/(L_g/k_g+R_c+L_s/k_s)\). This is a verified 1-D **conjugate/multiregion heat-transfer** unit, not an integrated 2-D OpenFOAM conjugate solver.
