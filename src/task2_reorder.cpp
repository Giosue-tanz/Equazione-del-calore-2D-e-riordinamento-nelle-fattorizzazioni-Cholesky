#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
using namespace std;

// Cartella di output condivisa con gli altri task
static const string OUTPUT_DIR = "output/";



struct Node {
    int id;  // identificatore naturale φ(i,j) = (i-1)*N + (j-1)
    int i;   // indice di riga   (1 ≤ i ≤ N)
    int j;   // indice di colonna (1 ≤ j ≤ N)
};

// semplicemente perche usare std::pair<int,int> sarebbe meno leggibile

// ============================================================
//  nested_dissection_geometrica
//  ===========================================================
//  Applica Divide et Impera partizionando geometricamente i nodi.
//  Complessità: O(N² log N) in tempo, profondità albero O(log N).
// ============================================================

// Oss:  porzione di griglia corrente da dividere, selettore per l'asse di taglio (0 = orizzontale/righe, 1 = verticale/colonne), lista globale condivisa (passata con & per poterla modificare senza crearne permutazione finale dei nodi che stiamo costruendo in post-ordine.

void nested_dissection_geometrica(const vector<Node>& nodi_iniziali, int axis_iniziale, vector<int>& ordering)
{
    // ── Caso base: nessun nodo ──────────────────────────────────────────
    if (nodi_iniziali.empty()) return;

    // ── Caso base: nodo singolo ─────────────────────────────────────────
    if (nodi_iniziali.size() == 1) {
        ordering.push_back(nodi_iniziali[0].id);
        return;
    }

    // ── Calcolo separatore geometrico ───────────────────────────────────
    int min_val = (axis_iniziale == 0) ? nodi_iniziali[0].i : nodi_iniziali[0].j; //op.ternario 
    int max_val = min_val;
    for (const auto& node : nodi_iniziali) {
        int val = (axis_iniziale == 0) ? node.i : node.j;
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
    }

    int mid_val = (min_val + max_val) / 2;

    vector<Node> V1, V2, VS;
    V1.reserve(nodi_iniziali.size());
    V2.reserve(nodi_iniziali.size()); // r evita riallocazioni multiple
    VS.reserve(nodi_iniziali.size());

    for (const auto& node : nodi_iniziali) { // nuova nomenclatura 
        int val = (axis_iniziale == 0) ? node.i : node.j; //op.ternario
        if      (val < mid_val) V1.push_back(node);
        else if (val > mid_val) V2.push_back(node);
        else                    VS.push_back(node);
    }

    // ── Chiamate ricorsive ──────────────────────────────────────────────
    // L'asse si alterna ad ogni livello di ricorsione.
    int next_axis = 1 - axis_iniziale;
    
    nested_dissection_geometrica(V1, next_axis, ordering);
    nested_dissection_geometrica(V2, next_axis, ordering);
    
    // ── Aggiunta separatore ─────────────────────────────────────────────
    // Il separatore VS viene aggiunto alla fine per essere posizionato in fondo
    for (const auto& node : VS) {
        ordering.push_back(node.id);
    }
}


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
    // Leggiamo solo id, i, j. Le coordinate x e y vengono ignorate per performance.
    while (coords_file >> n >> i >> j) {
        nodes.push_back({n, i, j});
        coords_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
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
