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
Il software adotta un'architettura modulare incentrata sul concetto di un risolutore generale denominato `HeatEquationSolver`. Di seguito vengono identificati i principali moduli che comporranno il sistema e i rispettivi compiti:

* **Modulo 1: Generazione della Griglia e Grafo di Adiacenza (C++)**
  * **Compito:** Riceve in input il parametro di discretizzazione $N$, genera i punti interni e calcola le relazioni di adiacenza orizzontale e verticale. Esporta i dati nei file di testo dedicati.
  
* **Modulo 2: Calcolo dell'Ordinamento Geometrico (C++)**
  * **Compito:** Implementa l'algoritmo di partizionamento ricorsivo dei nodi basato sulle coordinate spaziali, alternando la separazione lungo gli assi $x$ e $y$ (*nested dissection*) per generare un ordinamento ottimizzato.
  
* **Modulo 3: Costruzione del Sistema Lineare ed Esportazione (C++)**
  * **Compito:** Applica lo schema di discretizzazione dell'equazione differenziale, gestisce la sorgente di calore e incorpora le condizioni al contorno, assemblando la matrice sparsa e il vettore dei termini noti secondo l'ordinamento scelto.
  
* **Modulo 4: Risoluzione e Analisi Numerica (Python / Jupyter Notebook)**
  * **Compito:** Carica i file della matrice e del vettore dei termini noti, converte la struttura nel formato sparso CSC e calcola la fattorizzazione di Cholesky.

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
