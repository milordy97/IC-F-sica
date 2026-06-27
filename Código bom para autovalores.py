import numpy as np
import matplotlib.pyplot as plt
import unicodeit

gama_min = "γ"  # só para pegar o caracter
gama_mai = "Γ"



# parâmetross
N = 1000
rho = 0.5
Delta = 0.0
Gamma = 1.0
k = 1.0

### cadeia ##
L = N / rho
r = np.sort(np.random.uniform(0, L, N))

#  MATRIZ M
M = np.zeros((N, N), dtype=complex)
for j in range(N):
    for m in range(N):
        if j == m:
            M[j, m] = 1j* Delta - Gamma / 2.0
        else:
            dist = np.abs(r[j] - r[m])
            M[j, m] =  -(Gamma / 2.0) * np.exp(1j * k * dist)

### autovalores  ###
autovalores = np.linalg.eigvals(M)
autovalores_invertidos = -autovalores

###  real e imaginária ###
parte_real = np.real(autovalores_invertidos)  
parte_imag = np.imag(autovalores_invertidos)  

# sem log(0) ##
epsilon = 1e-15
parte_real_pos = np.abs(parte_real) + epsilon


plt.figure(figsize=(10, 6))

### Plotar  autovalores ##
plt.scatter(parte_imag, parte_real_pos, alpha=0.7, c='red', edgecolors='black', s=60, 
            label=f'Autovalores (N={N})')

# Linhas p/ referências
plt.axhline(y=epsilon, color='r', linestyle='--', alpha=0.5, label='$10^{-15}$')
plt.axhline(y=1.0, color='g', linestyle='--', alpha=0.5, label='$10^{0}$')
plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)


#  escala logarítmica no eixo y ###
plt.yscale('log')

# Labels e título
plt.xlabel('(\u03c9\u2099)', fontsize=14)
plt.ylabel('(γ\u2099)', fontsize=14)
plt.title(f' Autospectro de Autovalores de 1000 átomos', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.3)
plt.legend(loc='best')


#y_min = epsilon
#y_max = max(max(parte_real_pos)*1.1, 10)  # Garantir que 10⁰ (1.0) apareça bem
#plt.ylim(y_min, y_max)


plt.tight_layout()
plt.show()
