import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt


# physical parameters
A = 2.0
alpha = 1.0

# Spatial grid
L = 30
N = 3000
x_min = -L
x_max = L
x = np.linspace(x_min, x_max, N)
h = (x_max - x_min) / (N - 1)

# Breaking potentials
def V1(x, A, alpha, lamb):
    sech_sq = (1.0 / np.cosh(alpha * x))**2
    tanh_sq = np.tanh(alpha * x)**2
    tanh = np.tanh(alpha * x)
    
    W_sq = A**2 * tanh_sq + 2*A*lamb * tanh + lamb**2
    W_derv = A*alpha*sech_sq
    
    return W_sq - W_derv

def V2(x, A, alpha, lamb):
    sech_sq = (1.0 / np.cosh(alpha * x))**2
    tanh_sq = np.tanh(alpha * x)**2
    tanh = np.tanh(alpha * x)
    
    W_sq = A**2 * tanh_sq + 2*A*lamb * tanh + lamb**2
    W_derv = A*alpha*sech_sq
    
    return W_sq + W_derv

# Constants
diag_principal_cinetica = 2.0 / h**2
diag_fuera = -1.0 / h**2 * np.ones(N - 1)

# Perturbations array (depends of the size of A)
lambdas = np.linspace(0, 3.0, 15)

print(f"{'Lambda':<10} | {'E_0 (H1)':<12} | {'E_0 (H2)':<12}")
print("-" * 40)
E1_lista = []
E2_lista = []
# Susy rupture cycle
for lamb_val in lambdas:
    V1_array = V1(x, A, alpha, lamb_val)
    V2_array = V2(x, A, alpha, lamb_val)
    
    diag_ppal_1 = diag_principal_cinetica + V1_array
    diag_ppal_2 = diag_principal_cinetica + V2_array

    # Hamiltonian matrices as sparce tridiagonal
    H1 = diags([diag_fuera, diag_ppal_1, diag_fuera], [-1, 0, 1])
    H2 = diags([diag_fuera, diag_ppal_2, diag_fuera], [-1, 0, 1])
    
    # Diagonalization
    # Using 'SA' to get the smallest eigenvalue (ground state)
    E1, psi1 = eigsh(H1, k=2, which='SA')
    E2, psi2 = eigsh(H2, k=2, which='SA')
    E1_lista.append(E1[0])
    E2_lista.append(E2[0])
    # Ground state (index 0)
    print(f"{lamb_val:<10.4f} | {E1[0]:<12.6f} | {E2[0]:<12.6f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# First graph: Energy levels vs lambda (breaking parameter)
ax1.plot(lambdas, E1_lista, 'b-o', label=r'$E_0^{(1)}$ (Ground State H1)', linewidth=2)
ax1.plot(lambdas, E2_lista, 'r-s', label=r'$E_0^{(2)}$ (Ground State H2)', linewidth=2)
ax1.axvline(x=2.0, color='k', linestyle='--', label=r'Critical Threshold $\lambda = A$')

ax1.set_title('Spontaneous Supersymmetry Breaking', fontsize=14)
ax1.set_xlabel(r'Perturbation Parameter $\lambda$', fontsize=12)
ax1.set_ylabel('Energy', fontsize=12)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)

# Second graph: Potentials V1 and 2 for specific lambda value
lamb_plot = 3.0
V1_plot = V1(x, A, alpha, lamb_plot)
V2_plot = V2(x, A, alpha, lamb_plot)

ax2.plot(x, V1_plot, 'b-', label=r'$V_1(x)$', linewidth=2)
ax2.plot(x, V2_plot, 'r--', label=r'$V_2(x)$', linewidth=2)
ax2.set_xlim(-10, 5) 
ax2.set_ylim(-5, 20)
# 
ax2.set_title(r"Potentials' spatial collapse  ($\lambda = 3.0$)", fontsize=14)
ax2.set_xlabel('Position (x)', fontsize=12)
ax2.set_ylabel('Potential V(x)', fontsize=12)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()