# Equazione del calore 2D e riordinamento nelle fattorizzazioni Cholesky

*Progetto di Informatica con Laboratorio A.A. 2025/2026. Risoluzione dell'equazione del calore 2D e ordinamento nested dissection per fattorizzazione Cholesky.*

**Anno Accademico:** 2025/2026  
**Autore:** Giosuè Aiello  

---

## 1. Descrizione Generale e Obiettivi del Progetto
Il progetto si propone di risolvere numericamente il problema differenziale dell'equazione del calore 2D per determinare la temperatura di equilibrio su una piastra quadrata. Il dominio viene discretizzato mediante una griglia di punti interni ed esterni. 

L'obiettivo principale è analizzare l'efficacia del riordinamento dei nodi del grafo di adiacenza tramite un'euristica geometrica ricorsiva (*nested dissection*) rispetto all'ordinamento naturale, valutando l'impatto sulla successiva fattorizzazione di Cholesky del sistema lineare risultante.

---

## 2. Architettura del Software e Moduli Principali
Il software adotta un'architettura modulare incentrata sul concetto di un risolutore generale denominato `HES--HeatEquationSolver`. L'architettura prevede la divisione del problema in quattro moduli principali, di cui tre sviluppati in C++ e uno in Python, che comunicano in modo sequenziale scambiandosi file di testo.

**Modulo 1: Generazione della Griglia e Grafo di Adiacenza (`task1_grid.cpp`)**
Questo primo modulo funge da punto di ingresso. Riceve in input da riga di comando il parametro di discretizzazione $N$ e si occupa di generare la griglia di nodi fisici interni al dominio spaziale e il relativo grafo di adiacenza. Sfruttando liste di adiacenza implementate tramite `std::vector`, il modulo mappa la struttura e salva le informazioni in output sui file `coords.txt` e `connectivity.txt`. L'intera generazione richiede tempo e spazio lineare $\mathcal{O}(N^2)$ rispetto al numero totale dei nodi interni.

**Modulo 2: Calcolo dell'Ordinamento Geometrico (`task2_reorder.cpp`)**
È il nucleo algoritmico del software. Riceve in input il file delle coordinate `coords.txt` (ignorando appositamente il grafo di adiacenza per ottimizzare gli accessi in memoria) per calcolare una permutazione ottimizzata dei nodi. L'approccio si basa sulla *nested dissection geometrica*, un algoritmo divide-et-impera che partiziona il dominio separando i nodi rispetto al piano mediano, alternando gli assi $x$ e $y$. Utilizza strutture dati leggere (come vettori e struct ad hoc) e restituisce in output il file `ordering.txt`. L'algoritmo ha una complessità temporale attesa di $\mathcal{O}(N^2 \log N)$ e un'occupazione spaziale pari a $\mathcal{O}(N^2)$.

**Modulo 3: Costruzione del Sistema Lineare ed Esportazione (`task3_assemble.cpp`)**
Questo modulo applica lo schema di discretizzazione tramite stencil a 5 punti per trasformare il problema differenziale in un sistema algebrico. Riceve in input $N$, le coordinate, l'ordinamento (attivabile tramite il flag opzionale `-r`) e gestisce le eventuali condizioni di Dirichlet al contorno e la sorgente di calore. Sfruttando array di permutazione inversa, assembla la matrice sparsa simmetrica in una lista di triplette e calcola il vettore dei termini noti. Genera in output i file `A.txt` e `rhs.txt` operando nel rispetto della complessità asintotica ottimale $\mathcal{O}(N^2)$ sia in tempo che in spazio.

**Modulo 4: Risoluzione e Analisi Numerica (`solver.ipynb`)**
L'ultimo modulo è sviluppato interamente in Python tramite Jupyter Notebook e unifica la risoluzione matematica all'analisi sperimentale. Legge i file `A.txt` e `rhs.txt` per vari $N$, assembla la matrice nel formato compresso CSC di SciPy e calcola la soluzione esatta del sistema lineare servendosi della fattorizzazione di Cholesky. Oltre alla computazione vera e propria, questo ambiente genera automaticamente i grafici comparativi (come gli spy plot per la struttura della matrice) e la visualizzazione delle temperature. In termini di efficienza, l'ordinamento ottimizzato riduce drasticamente l'occupazione spaziale del fill-in da $\mathcal{O}(N^3)$ a $\mathcal{O}(N^2 \log N)$.

---

## 3. Interazioni tra i Moduli
I moduli interagiscono in modo sequenziale attraverso file di testo intermedi che fungono da interfacce di comunicazione:

1. Il **Modulo 1** produce i file `coords.txt` e `connectivity.txt`.
2. Il **Modulo 2** legge tali file e genera il file `ordering.txt` contenente la sequenza riordinata dei nodi.
3. Il **Modulo 3** utilizza sia le informazioni strutturali del grafo sia l'ordinamento per calcolare e produrre i file finali del sistema lineare, ossia `A.txt` e `rhs.txt`.
4. Il **Modulo 4** acquisisce `A.txt` e `rhs.txt` per eseguire la computazione finale e l'analisi sperimentale.

---

## 4. Strutture Dati Fondamentali
Nelle fasi iniziali dello sviluppo si prevede l'utilizzo delle seguenti strutture dati fondamentali:

* **Liste di adiacenza (`std::vector` in C++):** Utilizzate per memorizzare in modo efficiente i nodi e gli archi del grafo di discretizzazione del dominio.
* **Strutture dati di mappatura geometrica:** Contenitori di supporto utili a convertire in tempo lineare le coordinate spaziali o gli indici logici nei corrispondenti identificatori progressivi dei nodi.
* **Formato Sparso CSC (Compressed Sparse Column in SciPy/Python):** Struttura dati ottimizzata impiegata per la memorizzazione della matrice del sistema e per l'interfacciamento con la libreria di fattorizzazione CHOLMOD.
