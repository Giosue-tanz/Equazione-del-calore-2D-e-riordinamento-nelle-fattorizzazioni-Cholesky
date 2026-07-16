// ============================================================
//  task2_reorder.cpp
//  Modulo 2 — Nested Dissection GEOMETRICA
//
//  Calcola una permutazione dei nodi della griglia N×N tramite
//  l'algoritmo di nested dissection BASATO ESCLUSIVAMENTE sulle
//  coordinate spaziali (i, j) lette da coords.txt.
//
//  APPROCCIO GEOMETRICO (non su grafo):
//  ─────────────────────────────────────
//  Per griglie cartesiane regolari il separatore geometrico ottimale
//  coincide esattamente con il separatore topologico: la riga (o colonna)
//  mediana divide il grafo di adiacenza in due componenti bilanciate senza
//  archi tra loro, senza bisogno di esplorare la lista di adiacenza.
//
//  Vantaggi:
//    1. Non richiede la costruzione/esplorazione del grafo  → O(N²logN)
//    2. Accesso sequenziale alle coordinate                 → cache-friendly
//    3. Implementazione O(N²logN) deterministica e ottimale per griglie
//    4. Separatore garantito bilanciato (N/2 nodi per lato)
//
//  Algoritmo (iterativo con stack esplicito):
//    • Inizia con tutti i nodi, asse di taglio = X (i)
//    • Calcola mid = (min_val + max_val) / 2 sull'asse corrente
//    • Partiziona: V1 (coord < mid), V2 (coord > mid), VS (coord = mid)
//    • Push su stack: VS (separatore), V2, V1  [ordine LIFO]
//    • Ordine finale: V1 → V2 → VS (separatore in fondo)
//    • Alterna asse X/Y ad ogni livello
//
//  Input:  output/coords.txt   (prodotto da task1_grid)
//  Output: output/ordering.txt (una riga: "new_id old_id")
//
//  Uso: ./task2_reorder
// ============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <stack>    // stack esplicito in heap (evita stack overflow su griglie grandi)

using namespace std;

// Cartella di output condivisa con gli altri task
static const string OUTPUT_DIR = "output/";

// ============================================================
//  Struttura nodo: mantiene sia gli indici logici (i, j) che
//  le coordinate fisiche (x, y). L'algoritmo geometrico usa
//  gli indici interi (i, j) per il taglio mediano, che sono
//  più robusti di confronti in floating-point.
// ============================================================
struct Node {
    int id;  // identificatore naturale φ(i,j) = (i-1)*N + (j-1)
    int i;   // indice di riga   (1 ≤ i ≤ N)
    int j;   // indice di colonna (1 ≤ j ≤ N)
    double x; // coordinata fisica x = i/(N+1)
    double y; // coordinata fisica y = j/(N+1)
};

// ============================================================
//  Item sullo stack esplicito per la nested dissection iterativa.
//  Separare il flag 'separatore' consente di ritardare la scrittura
//  del separatore VS fino a che V1 e V2 non sono completati.
// ============================================================
struct StackItem {
    vector<Node> S;       // sottoinsieme di nodi da processare
    int axis;             // asse di taglio: 0 = X (colonne di i), 1 = Y (colonne di j)
    bool is_separator;    // true: questo item va aggiunto direttamente all'ordinamento
};

// ============================================================
//  nested_dissection_geometrica
//  ────────────────────────────
//  Versione ITERATIVA con stack esplicito allocato sull'heap.
//  La versione ricorsiva sarebbe più leggibile ma causerebbe
//  stack overflow per N >> 1000 (profondità O(log N) con sotto-
//  problemi di dimensione O(N²) porta a frame molto grandi).
//
//  Complessità: O(N² log N) in tempo, O(N²) in spazio.
// ============================================================
void nested_dissection_geometrica(const vector<Node>& nodi_iniziali,
                                   int axis_iniziale,
                                   vector<int>& ordering)
{
    stack<StackItem> stk;
    stk.push({nodi_iniziali, axis_iniziale, false});

    while (!stk.empty()) {
        StackItem item = move(stk.top());
        stk.pop();

        // ── Caso 1: separatore già calcolato → aggiungilo all'ordinamento ──
        if (item.is_separator) {
            for (const auto& node : item.S)
                ordering.push_back(node.id);
            continue;
        }

        // ── Caso base: nessun nodo ──────────────────────────────────────────
        if (item.S.empty()) continue;

        // ── Caso base: nodo singolo ─────────────────────────────────────────
        if (item.S.size() == 1) {
            ordering.push_back(item.S[0].id);
            continue;
        }

        // ── Calcolo separatore geometrico ───────────────────────────────────
        // Usiamo gli indici interi (i o j) per robustezza numerica.
        // Il taglio mediano garantisce bilanciamento esatto per griglie regolari.
        int min_val = (item.axis == 0) ? item.S[0].i : item.S[0].j;
        int max_val = min_val;
        for (const auto& node : item.S) {
            int val = (item.axis == 0) ? node.i : node.j;
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
        }

        // Mediana intera: separa [min_val, mid-1] da [mid+1, max_val]
        // VS = { nodi con coordinata == mid }  → separatore geometrico
        int mid_val = (min_val + max_val) / 2;

        vector<Node> V1, V2, VS;
        V1.reserve(item.S.size());
        V2.reserve(item.S.size());
        VS.reserve(item.S.size());

        for (const auto& node : item.S) {
            int val = (item.axis == 0) ? node.i : node.j;
            if      (val < mid_val) V1.push_back(node);
            else if (val > mid_val) V2.push_back(node);
            else                    VS.push_back(node);
        }

        // ── Push in ordine LIFO ─────────────────────────────────────────────
        // Ordine desiderato:   V1 → V2 → VS
        // Quindi push inverso: VS (fondo), V2, V1 (cima)
        // L'asse si alterna ad ogni livello di ricorsione.
        int next_axis = 1 - item.axis;
        stk.push({VS, item.axis, true});   // separatore: flag = true
        stk.push({V2, next_axis, false});
        stk.push({V1, next_axis, false});
    }
}

// ============================================================
//  main
// ============================================================
int main() {

    // ── 1. Lettura di coords.txt ────────────────────────────────────────────
    // Il Modulo 2 usa SOLO le coordinate geometriche dei nodi.
    // Il file connectivity.txt (grafo di adiacenza) NON viene letto:
    // per griglie regolari la struttura topologica è interamente codificata
    // dagli indici (i, j) e il separatore geometrico è anche separatore del grafo.
    ifstream coords_file(OUTPUT_DIR + "coords.txt");
    if (!coords_file.is_open()) {
        cerr << "Errore: impossibile aprire " << OUTPUT_DIR << "coords.txt" << endl;
        cerr << "Assicurarsi di aver eseguito task1_grid prima." << endl;
        return 1;
    }

    vector<Node> nodes;
    int n, i, j;
    double x, y;
    while (coords_file >> n >> i >> j >> x >> y) {
        nodes.push_back({n, i, j, x, y});
    }
    coords_file.close();

    if (nodes.empty()) {
        cerr << "Errore: nessun nodo trovato in " << OUTPUT_DIR << "coords.txt" << endl;
        return 1;
    }

    // ── 2. Nested Dissection Geometrica ────────────────────────────────────
    vector<int> ordering;
    ordering.reserve(nodes.size());

    // Primo taglio lungo X (asse 0 = indice di riga i)
    nested_dissection_geometrica(nodes, 0, ordering);

    // Verifica invariante: ordering deve essere una permutazione valida
    if ((int)ordering.size() != (int)nodes.size()) {
        cerr << "Errore critico: ordering ha " << ordering.size()
             << " elementi, attesi " << nodes.size() << endl;
        return 1;
    }

    // ── 3. Salvataggio dell'ordinamento ────────────────────────────────────
    // Formato: "new_id old_id"
    //   new_id  = posizione nel nuovo ordinamento (0-based)
    //   old_id  = identificatore naturale del nodo φ(i,j)
    ofstream out_file(OUTPUT_DIR + "ordering.txt");
    if (!out_file.is_open()) {
        cerr << "Errore: impossibile creare " << OUTPUT_DIR << "ordering.txt" << endl;
        return 1;
    }

    for (size_t k = 0; k < ordering.size(); ++k) {
        out_file << k << " " << ordering[k] << "\n";
    }
    out_file.close();

    // ── 4. Report ──────────────────────────────────────────────────────────
    cout << "Nested Dissection Geometrica completata." << endl;
    cout << "  Nodi elaborati : " << nodes.size() << endl;
    cout << "  File generato  : " << OUTPUT_DIR << "ordering.txt" << endl;

    return 0;
}
