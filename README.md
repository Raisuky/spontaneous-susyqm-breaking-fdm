# Numerical study of spontaneous SUSY breaking

## Overview
This repository contains a comp. physics project which demonstrates the spontaneous breaking of Supersymmetry (SUSY) in Quantum Mechanics. This was done by applying the Finite Difference Method (FDM), where we solved the stationary Schrödingers eq. using a modified asymmetric Pöschl-teller potential (Rosen-Morse II potential) and analizing the energy spectra of the partner hamiltonians.

## Theoretical background
In unbroken SUSY QM, the ground state of enegery of the first hamiltonian is zero $E_0^{(1)}$, whereas the spectra of the partner Hamiltonian $H_1$ and $H_2$ are completely degenerate.
To induce the spontaneous breaking we introduce an assymetric parameter $\lambda$ to the superpotential:
$$ W(x) = A \tanh(\alpha x) + \lambda$$

The theoretical constrains dictate that when $\lambda$ crosses the threshold $\lambda \ge A$ the groundstate is no longer normalizable, which forces the $E_0^{(1)} > 0$, breaking the suppersimetry.
## Methodology
- **Numerical Method:** Matrix diagonalization via FDM using a discrete spatial grid.
- **Tools:** Python, numpy, scipy.sparse, and matplotlib.
- **Simulation Parameters:** The spatial domain is truncated using Dirichlet boundary conditions at $x=\pm30$ with a grid resolution of $N = 3000$ points to minimize the finite-size effects and truncation errors.

## Key results
