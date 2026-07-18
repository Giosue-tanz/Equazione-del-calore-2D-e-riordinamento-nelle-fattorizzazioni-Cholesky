#!/usr/bin/env python
# coding: utf-8

# # Analisi Numerica del Problema di Poisson 2D: Ordinamento Naturale vs Nested Dissection
# 
# ## 1. Introduzione Teorica
# 
# In questa relazione tecnica si affronta la risoluzione numerica dell'Equazione del Calore stazionaria bidimensionale (Equazione di Poisson) su un dominio quadrato $\Omega = [0,1] \times [0,1]$:
# $$
# - \Delta u(x,y) = f(x,y) \quad \text{in } \Omega
# $$
# 
# soggetta a condizioni al contorno di Dirichlet $u(x,y) = u_0(x,y)$ su $\partial \Omega$. 
# Nello specifico, consideriamo assenza di forzante interna ($f(x,y) = 0$, equazione di Laplace) e una forzante al bordo armonica definita come $u_0(x,y) = \sin(\pi x)$. Come implementato nel modulo `task3_assemble.cpp`, questa condizione al contorno viene applicata spazialmente su tutti i bordi della griglia discretizzata, manifestando la massima energia alle coordinate $y=0$ e $y=1$ lungo l'asse $x$, annullandosi ai bordi laterali $x=0$ e $x=1$.
# 
# Per la discretizzazione spaziale si adotta il metodo delle **Differenze Finite (FDM)** con passo spaziale uniforme $h = \frac{1}{N+1}$. Questo produce un sistema di equazioni algebriche lineari del tipo:
# $$
# A x = b
# $$
# dove $A \in \mathbb{R}^{N^2 \times N^2}$ è una matrice **sparsa, simmetrica e definita positiva**, e $x$ rappresenta i valori discretizzati del campo $u$ nei punti interni della griglia. L'operatore differenziale (Laplaciano) è approssimato tramite stencil a 5 punti, per cui la matrice $A$ presenta una struttura a banda ben determinata.
# 
# Il sistema viene risolto tramite metodo diretto, basato sulla **Fattorizzazione di Cholesky**:
# $$
# A = L L^T
# $$
# con $L$ matrice triangolare inferiore. 
# L'obiettivo dell'analisi è misurare e confrontare il **fill-in** (perdita di sparsità) della matrice $L$ generata applicando l'ordinamento spaziale "Naturale" rispetto all'algoritmo algrebrico di "**Nested Dissection**", che mira a minimizzare il riempimento di banda.
# 

# In[1]:


import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import subprocess
import time
from pathlib import Path
from sksparse.cholmod import cholesky
from matplotlib import cm

# Disabilitazione di warning inutili per plot
import warnings
warnings.filterwarnings('ignore')

ROOT = Path.cwd()
if ROOT.name == 'notebook':
    ROOT = ROOT.parent

OUTPUT = ROOT / "output"
OUTPUT.mkdir(exist_ok=True)


# ## 2. Metodologia di Risoluzione
# 
# Implementiamo le routine per il parsing dei file contenenti la struttura della matrice sparsa e del vettore dei termini noti, generati dalla porzione C++ del progetto (Task 3).
# Sfruttiamo l'efficiente formato `CSC` (Compressed Sparse Column) che si adatta nativamente al solutore CHOLMOD.
# 

# In[2]:


def read_A(filename):
    """Legge le triplette (i, j, val) e costruisce la matrice in formato CSC."""
    data = np.loadtxt(filename)
    if data.size == 0:
        return sp.csc_matrix((0,0))
    rows = data[:, 0].astype(int)
    cols = data[:, 1].astype(int)
    vals = data[:, 2]
    # L'ordinamento è 0-based
    N_tot = max(rows.max(), cols.max()) + 1
    A = sp.coo_matrix((vals, (rows, cols)), shape=(N_tot, N_tot))
    return A.tocsc()

def read_b(filename):
    """Legge il vettore del termine noto."""
    return np.loadtxt(filename)

def run_pipeline(N, use_reorder=False, bordo=3):
    """
    Wrapper per l'esecuzione in background della pipeline C++:
    1. Genera griglia ed ordinamento (task1)
    2. Calcola connettività e matrice A (task3)
    """
    cmd1 = [str(ROOT / 'task1'), str(N)]
    subprocess.run(cmd1, cwd=str(ROOT), check=True, capture_output=True)

    if use_reorder:
        cmd2 = [str(ROOT / 'task2')]
        subprocess.run(cmd2, cwd=str(ROOT), check=True, capture_output=True)

    cmd3 = [str(ROOT / 'task3'), str(N)]
    if use_reorder:
        cmd3.append('-r')
    subprocess.run(cmd3, input=str(bordo).encode(), cwd=str(ROOT), check=True, capture_output=True)

def my_cholesky(A):
    """
    Esegue la fattorizzazione di Cholesky sfruttando scikit-sparse.
    CHOLMOD viene forzato ad utilizzare l'ordinamento "natural", 
    ovvero l'ordine imposto dalla matrice passata (Naturale o Nested Dissection che sia).
    Ritorna il fattore triangolare inferiore L in formato CSC.
    """
    factor = cholesky(A, order="natural")
    if isinstance(factor, tuple):
        return factor[0].tocsc() if hasattr(factor[0], "tocsc") else factor[0]
    return factor.L()


# ## 3. Analisi della Struttura Matriciale per $N=4$
# 
# Iniziamo la validazione strutturale per un dominio estremamente ridotto: una griglia interna di $4 	imes 4$ nodi, corrispondente a un sistema di $N_{incognite} = 16$. 
# La matrice coefficienti risulterà quindi di dimensioni $16 	imes 16$.
# Risolviamo inizialmente il sistema con la configurazione base (Ordinamento Naturale, Condizione al Bordo 3).
# 

# In[3]:


# Generazione griglia e assemblaggio matrice (N=4) con Ordinamento Naturale
N_test = 4
run_pipeline(N_test, use_reorder=False, bordo=3)

# Lettura dai file generati
A_n4 = read_A(OUTPUT / "A.txt")
b_n4 = read_b(OUTPUT / "rhs.txt")

# Calcolo del fattore di Cholesky di -A
# (Il Laplaciano discreto -A è Definito Positivo)
neg_A_n4 = (-A_n4).tocsc()
L_n4 = my_cholesky(neg_A_n4)

# Risoluzione del sistema L L^T x = b
# L y = -b  =>  y = L \ -b
y_n4 = spla.spsolve_triangular(L_n4, -b_n4, lower=True)
# L^T x = y =>  x = L^T \ y
u_n4 = spla.spsolve_triangular(L_n4.T.tocsr(), y_n4, lower=False)

# Visualizzazione (Spy Plot)
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Usiamo marker netti (quadrati) proporzionati per mostrare fedelmente gli elementi della matrice
ax[0].spy(A_n4, markersize=15, marker='s', color='navy')
ax[0].set_title(r"Matrice $A$ ($16 	imes 16$)")
ax[0].set_xticks(np.arange(0, 16, 4))
ax[0].set_yticks(np.arange(0, 16, 4))
ax[0].grid(True, linestyle='--', alpha=0.5)

ax[1].spy(L_n4, markersize=15, marker='s', color='darkred')
ax[1].set_title(r"Fattore $L$ ($16 	imes 16$)")
ax[1].set_xticks(np.arange(0, 16, 4))
ax[1].set_yticks(np.arange(0, 16, 4))
ax[1].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig(OUTPUT / "spy_N4.png", dpi=150)
plt.show()

print(f"Norma del Residuo ||A x - b|| per N={N_test}: {np.linalg.norm(A_n4.dot(u_n4) - b_n4):.2e}")


# ### Osservazioni sullo Spy Plot (Ordinamento Naturale)
# Osservando la **Matrice $A$** sulla sinistra, notiamo il classico pattern a blocchi tridiagonale indotto dallo stencil a 5 punti con ordinamento lessicografico. Sono ben distinguibili la diagonale principale e le sottomatrici di accoppiamento.
# Nel calcolo del **Fattore $L$** (destra), è visibile il fenomeno del **Fill-In**: lo spazio vuoto tra la diagonale e la "sotto-banda" più esterna (distante $N=4$ posizioni) viene interamente saturato da elementi non nulli, dovuti alle operazioni algebriche di eliminazione. Per griglie di grandi dimensioni ($N$ elevato), la banda diventerà enormemente densa, rendendo la memorizzazione in RAM insostenibile $O(N^3)$.
# 

# ## 4. Ricostruzione del Campo Termico (Soluzione Bidimensionale)
# 
# Il vettore soluzione $x$ (o `u_n4`) corrisponde al campo scalare mappato vettorialmente in forma un-dimensionale `(16 x 1)`. 
# Provvediamo al reshaping matriciale $N 	imes N$ includendo geometricamente i punti di frontiera, in modo da poter visualizzare la temperatura in tutto il dominio spaziale $\Omega$.
# 

# In[4]:


# Rigeneriamo la griglia per N=32 al fine di ottenere un grafico più morbido
N_vis = 32
run_pipeline(N_vis, use_reorder=False, bordo=3)
A_vis = read_A(OUTPUT / "A.txt")
b_vis = read_b(OUTPUT / "rhs.txt")

L_vis = my_cholesky((-A_vis).tocsc())
y_vis = spla.spsolve_triangular(L_vis, -b_vis, lower=True)
u_vis = spla.spsolve_triangular(L_vis.T.tocsr(), y_vis, lower=False)

# Reshape del vettore soluzione in una griglia 2D (N x N)
U_inner = u_vis.reshape((N_vis, N_vis))

# Estensione del dominio con i bordi (N+2 x N+2)
U_full = np.zeros((N_vis + 2, N_vis + 2))
U_full[1:-1, 1:-1] = U_inner

# Applicazione della boundary condition u0 = sin(pi*x)
# In accordo con get_u0, la forzante dipende spazialmente solo dalla coordinata X 
# ma è geometricamente imposta su tutte le posizioni fisiche di frontiera.
x_coords = np.linspace(0, 1, N_vis + 2)
boundary_vals = np.sin(np.pi * x_coords)

U_full[0, :] = boundary_vals   # Bordo Y=0 (bottom)
U_full[-1, :] = boundary_vals  # Bordo Y=1 (top)
U_full[:, 0] = 0.0             # Bordo X=0: sin(0) = 0
U_full[:, -1] = 0.0            # Bordo X=1: sin(pi) = 0

# Generazione del plot 3D del Campo Termico
X, Y = np.meshgrid(x_coords, np.linspace(0, 1, N_vis + 2))

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, U_full, cmap='plasma', edgecolor='none')
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Temperatura u(x,y)')

ax.set_title(r"Campo Termico 3D $u(x,y)$ per $N=32$")
ax.set_xlabel(r"Asse $x$")
ax.set_ylabel(r"Asse $y$")
ax.set_zlabel(r"Temperatura $u$")

plt.savefig(OUTPUT / "heatmap_N32.png", dpi=150)
plt.show()


# L'analisi della mappa di calore mostra rigorosamente la coerenza con le equazioni di Dirichlet:
# 1. L'energia immessa nel sistema si localizza in $Y=0$ e $Y=1$ formando i due lobi ad alta temperatura, derivati dall'applicazione della condizione analitica $u_0 = \sin(\pi x)$.
# 2. I gradienti termici convergono verso le pareti omogenee $X=0$ e $X=1$ e sfumano verso il centro geometrico del dominio in maniera totalmente simmetrica.
# 

# ## 5. Impatto del Riordinamento (Nested Dissection)
# 
# Utilizziamo ora la pipeline su una griglia di dimensioni $32 	imes 32$ ($N_{incognite}=1024$) attivando l'algoritmo di **Nested Dissection** per valutarne le proprietà. Tale riordinamento divide il grafo della matrice tramite l'identificazione iterativa di grafi separatori (Vertex Separators), riorganizzando la struttura della matrice in macro-blocchi gerarchici, frammentando la banda originaria per limitare la densità di riempimento.
# 

# In[5]:


# Generazione con Nested Dissection
run_pipeline(N_vis, use_reorder=True, bordo=3)
A_nd = read_A(OUTPUT / "A.txt")
L_nd = my_cholesky((-A_nd).tocsc())

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Plot L Naturale (calcolato nel blocco precedente)
axs[0].spy(L_vis, markersize=0.5, color='darkred')
axs[0].set_title(f"L (Ordinamento Naturale)\nnnz = {L_vis.nnz:,}")

# Plot L Nested Dissection
axs[1].spy(L_nd, markersize=0.5, color='darkgreen')
axs[1].set_title(f"L (Nested Dissection)\nnnz = {L_nd.nnz:,}")

plt.tight_layout()
plt.savefig(OUTPUT / "spy_N32.png", dpi=150)
plt.show()

print(f"Rapporto di Fill-in L_Nat / L_ND = {L_vis.nnz / L_nd.nnz:.2f}x")


# ## 6. Benchmark Prestazionale Asintotico
# 
# Passiamo ora a testare empiricamente la scalabilità algoritmica di entrambe le strategie di ordinamento per una sequenza di discretizzazioni via via più fitte: $N \in \{32, 64, 128, 256, 512, 1024\}$.
# **Nota Accademica**: L'ordinamento *Naturale* per $N=1024$ richiederà oltre $\approx 13\text{ GB}$ di RAM per ospitare oltre $1$ miliardo di zeri diventati non nulli, generando una potenziale rottura per *Out Of Memory*. Per questo motivo l'ordinamento naturale verrà interrotto a $N=512$, mentre il Nested Dissection procederà agilmente alla massima discretizzazione in virtù dell'efficienza della banda frammentata.
# 

# In[6]:


Ns = [32, 64, 128, 256, 512, 1024]

res_nat = {"nnz_A": [], "nnz_L": [], "t_chol": []}
res_nd  = {"nnz_A": [], "nnz_L": [], "t_chol": []}

for N in Ns:
    print(f"\n{'='*60}")
    print(f"  Analisi per N = {N}  (Sistema: {N*N} x {N*N})")
    print(f"{'='*60}")

    for label, use_reorder, dataset in [("Naturale", False, res_nat), 
                                        ("Nested Dissection", True, res_nd)]:

        # Evitiamo OOM sul Naturale a 1024
        if label == "Naturale" and N >= 1024:
            print(f"  [{label:<17}]  -- Skipped per restrizioni di Memoria RAM --")
            continue

        run_pipeline(N, use_reorder=use_reorder)
        A_cur = read_A(OUTPUT / "A.txt")
        neg_A = (-A_cur).tocsc()

        t0 = time.perf_counter()
        L_cur = my_cholesky(neg_A)
        t_chol = time.perf_counter() - t0

        dataset["nnz_A"].append(A_cur.nnz)
        dataset["nnz_L"].append(L_cur.nnz)
        dataset["t_chol"].append(t_chol)

        print(f"  [{label:<17}]  nnz(A) = {A_cur.nnz:>9,}  |  "
              f"nnz(L) = {L_cur.nnz:>11,}  |  "
              f"T_chol = {t_chol:>7.3f} s")


# In[7]:


fig, ax = plt.subplots(1, 2, figsize=(14, 5))

# Estrazione array per calcoli
N_nat = np.array(Ns[:len(res_nat['nnz_L'])])
N_nd = np.array(Ns[:len(res_nd['nnz_L'])])
dof_nat = N_nat**2
dof_nd = N_nd**2

# Plot Fill-in
ax[0].loglog(dof_nat, res_nat["nnz_L"], 's-', color='darkred', label="Naturale", linewidth=2)
ax[0].loglog(dof_nd, res_nd["nnz_L"], 'o-', color='darkgreen', label="Nested Dissection", linewidth=2)
# Curve Teoriche
ax[0].loglog(dof_nd, dof_nd**(3/2), '--', color='gray', label=r"Teorica $O(N_{dof}^{3/2})$")
ax[0].loglog(dof_nd, dof_nd * np.log(dof_nd), ':', color='black', label=r"Teorica $O(N_{dof} \log N_{dof})$")

ax[0].set_title(r"Analisi Asintotica del Fill-in: $nnz(L)$")
ax[0].set_xlabel("Gradi di Libertà ($N^2$)")
ax[0].set_ylabel("Numero elementi non-zero in L")
ax[0].legend()
ax[0].grid(True, which="both", ls="--", alpha=0.5)

# Plot Tempi
ax[1].loglog(dof_nat, res_nat["t_chol"], 's-', color='darkred', label="Naturale", linewidth=2)
ax[1].loglog(dof_nd, res_nd["t_chol"], 'o-', color='darkgreen', label="Nested Dissection", linewidth=2)

ax[1].set_title("Costo Computazionale: Tempo Fattorizzazione Cholesky")
ax[1].set_xlabel("Gradi di Libertà ($N^2$)")
ax[1].set_ylabel("Tempo (s)")
ax[1].legend()
ax[1].grid(True, which="both", ls="--", alpha=0.5)

plt.tight_layout()
plt.savefig(OUTPUT / "benchmark.png", dpi=150)
plt.show()


# ## 7. Discussione dei Risultati e Conclusioni
# 
# L'evidenza sperimentale corrobora perfettamente la teoria analitica delle matrici sparse per FDM 2D:
# 
# 1. **Andamento dell'Ordinamento Naturale**: Il numero di zeri che subisce *fill-in* cresce come $\sim O(N_{dof}^{1.5})$ ovvero $O(N^3)$, a causa della semi-banda di dimensione $N$. Di conseguenza, i tempi di esecuzione degradano molto rapidamente, portando a limitazioni ferree per griglie che superano i $512 	imes 512$ nodi.
# 2. **Andamento di Nested Dissection**: La divisione ricorsiva del dominio produce una matrice di permutazione che limita drasticamente l'operatore di eliminazione, portando la densità dei non-zero nel fattore $L$ a scalare approssimativamente in $\sim O(N_{dof} \log N_{dof})$. Ciò si traduce in una curva di crescita quasi lineare nel plot logaritmico rispetto all'estrema severità strutturale dell'ordinamento lessicografico.
# 3. Il fattore **prestazionale (CPU Time)** segue linearmente i vantaggi strutturali misurati, dimostrando una riduzione nel tempo computazionale che sorpassa ampiamente il fattore di guadagno $100	imes$ già su reti moderatamente complesse ($512^2$), permettendo di gestire senza fatica l'imponente sistema a $1024^2$ gradi di libertà (oltre 1 milione di incognite) in pochi secondi.
# 
