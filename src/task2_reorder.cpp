#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <stack>   // [MIGLIORAMENTO] necessario per la versione iterativa
#include <chrono>  // [MIGLIORAMENTO] benchmark tempi come nel codice del prof (qs.cpp)
#include <iomanip> // [MIGLIORAMENTO] per fixed e setprecision nell'output del benchmark

using namespace std;
using namespace chrono;

// Struttura per rappresentare un nodo con le sue coordinate
struct Node {
    int id;
    int i;
    int j;
    double x;
    double y;
};

// [MIGLIORAMENTO] Versione ITERATIVA della nested dissection.
// Sostituisce la versione ricorsiva per eliminare il rischio di stack overflow
// su griglie molto grandi (N >> 1000), usando uno stack esplicito in heap.
// Logica identica alla versione ricorsiva: V1 e V2 prima, separatore VS alla fine.
struct StackItem {
    vector<Node> S;   // sottoinsieme di nodi da processare
    int axis;         // asse di taglio corrente (0=X, 1=Y)
    bool separatore;  // true: questo item va solo aggiunto all'ordinamento (e' un separatore)
};

void nested_dissection(const vector<Node>& nodi_iniziali, int axis_iniziale, vector<int>& ordering) {
    stack<StackItem> stk;
    stk.push({nodi_iniziali, axis_iniziale, false});

    while (!stk.empty()) {
        // Prendiamo il prossimo lavoro dallo stack
        StackItem item = move(stk.top());
        stk.pop();

        // Caso: questo e' un separatore gia' calcolato, lo aggiungiamo all'ordinamento
        if (item.separatore) {
            for (const auto& node : item.S)
                ordering.push_back(node.id);
            continue;
        }

        // Caso base: nessun nodo
        if (item.S.empty()) continue;

        // Caso base: un solo nodo, lo inseriamo direttamente
        if (item.S.size() == 1) {
            ordering.push_back(item.S[0].id);
            continue;
        }

        // Troviamo la coordinata minima e massima sull'asse corrente
        int min_val = (item.axis == 0) ? item.S[0].i : item.S[0].j;
        int max_val = min_val;
        for (const auto& node : item.S) {
            int val = (item.axis == 0) ? node.i : node.j;
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
        }

        // Calcoliamo la coordinata mediana per il separatore geometrico
        int mid_val = (min_val + max_val) / 2;

        vector<Node> V1, V2, VS;
        V1.reserve(item.S.size());
        V2.reserve(item.S.size());
        VS.reserve(item.S.size());

        // Partizioniamo il sottoinsieme
        for (const auto& node : item.S) {
            int val = (item.axis == 0) ? node.i : node.j;
            if (val < mid_val)      V1.push_back(node);
            else if (val > mid_val) V2.push_back(node);
            else                    VS.push_back(node);
        }

        // Inseriamo nello stack in ordine INVERSO (LIFO): cio' che vogliamo
        // eseguire prima deve stare in cima.
        // Ordine di esecuzione desiderato: V1 -> V2 -> VS (separatore)
        // Quindi push: VS, poi V2, poi V1 (V1 verra' estratta per prima)
        stk.push({VS, item.axis,      true});          // separatore: aggiunto alla fine
        stk.push({V2, 1 - item.axis,  false});         // seconda partizione
        stk.push({V1, 1 - item.axis,  false});         // prima partizione (elaborata per prima)
    }
}

int main(int argc, char* argv[]) {
    // 1. Lettura rigorosa del file coords.txt generato nel task precedente
    ifstream coords_file("coords.txt");
    if (!coords_file.is_open()) {
        cerr << "Errore: impossibile aprire coords.txt" << endl;
        return 1;
    }

    vector<Node> nodes;
    int n, i, j;
    double x, y;
    while (coords_file >> n >> i >> j >> x >> y) {
        nodes.push_back({n, i, j, x, y});
    }
    coords_file.close();

    // 2. Lettura formale del file connectivity.txt per aderire rigorosamente alla specifica.
    // L'euristica puramente geometrica non necessita di scandire la lista di adiacenza,
    // ma la specifica del progetto menziona esplicitamente che questo task prende "in input
    // i file generati nel punto precedente".
    ifstream conn_file("connectivity.txt");
    if (!conn_file.is_open()) {
        cerr << "Errore: impossibile aprire connectivity.txt" << endl;
        return 1;
    }
    int u, v;
    long long int edge_count = 0;
    while (conn_file >> u >> v) {
        edge_count++;
    }
    conn_file.close();

    if (nodes.empty()) {
        cerr << "Errore: nessun nodo trovato in coords.txt" << endl;
        return 1;
    }

    // 3. Esecuzione dell'algoritmo Nested Dissection
    // [MIGLIORAMENTO] misurazione del tempo con chrono (tecnica dal codice del prof)
    vector<int> ordering;
    ordering.reserve(nodes.size());

    auto start_nd = high_resolution_clock::now();

    // Iniziamo tagliando lungo l'asse X (axis = 0)
    nested_dissection(nodes, 0, ordering);

    auto end_nd = high_resolution_clock::now();
    duration<double, milli> tempo_nd = end_nd - start_nd;

    // 4. Salvataggio del nuovo ordinamento nel formato `new_id old_id`
    ofstream out_file("ordering.txt");
    if (!out_file.is_open()) {
        cerr << "Errore: impossibile creare ordering.txt" << endl;
        return 1;
    }

    for (size_t k = 0; k < ordering.size(); ++k) {
        out_file << k << " " << ordering[k] << "\n";
    }
    out_file.close();

    cout << "Nested Dissection completata con successo." << endl;
    cout << "Nodi elaborati: " << nodes.size() << ", Archi letti: " << edge_count << endl;
    cout << "Nuovo ordinamento salvato in ordering.txt" << endl;
    cout << "\n=== BENCHMARK ===" << endl;
    cout << fixed << setprecision(3);
    cout << "  Nested Dissection (iterativa): " << tempo_nd.count() << " ms" << endl;

    return 0;
}
