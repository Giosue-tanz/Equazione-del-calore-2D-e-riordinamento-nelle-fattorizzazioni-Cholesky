#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

// Struttura per rappresentare un nodo con le sue coordinate
struct Node {
    int id;
    int i;
    int j;
    double x;
    double y;
};

// Funzione ricorsiva per la nested dissection
// S: sottoinsieme dei nodi da processare
// axis: 0 per tagliare lungo l'asse X (usando i), 1 per tagliare lungo l'asse Y (usando j)
// ordering: lista globale in cui accumulare l'ordinamento finale (old_id)
void nested_dissection(const vector<Node>& S, int axis, vector<int>& ordering) {
    if (S.empty()) {
        return;
    }
    if (S.size() == 1) {
        ordering.push_back(S[0].id);
        return;
    }

    // Troviamo la coordinata minima e massima sull'asse corrente
    int min_val = (axis == 0) ? S[0].i : S[0].j;
    int max_val = min_val;

    for (const auto& node : S) {
        int val = (axis == 0) ? node.i : node.j;
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
    }

    // Calcoliamo la coordinata mediana per il separatore geometrico
    int mid_val = (min_val + max_val) / 2;

    vector<Node> V1, V2, VS;
    // Preallochiamo memoria per evitare riallocazioni continue
    V1.reserve(S.size());
    V2.reserve(S.size());
    VS.reserve(S.size());

    // Partizioniamo il sottoinsieme S
    for (const auto& node : S) {
        int val = (axis == 0) ? node.i : node.j;
        if (val < mid_val) {
            V1.push_back(node);
        } else if (val > mid_val) {
            V2.push_back(node);
        } else {
            VS.push_back(node);
        }
    }

    // Chiamata ricorsiva sulle due partizioni invertendo l'asse di taglio
    nested_dissection(V1, 1 - axis, ordering);
    nested_dissection(V2, 1 - axis, ordering);

    // I nodi del separatore (VS) vengono inseriti alla fine dell'ordinamento
    // per posticipare i loro archi (e limitare il fill-in) durante la fattorizzazione
    for (const auto& node : VS) {
        ordering.push_back(node.id);
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
    vector<int> ordering;
    ordering.reserve(nodes.size());

    // Iniziamo tagliando lungo l'asse X (axis = 0)
    nested_dissection(nodes, 0, ordering);

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
    cout << "Nodi elaborati: " << nodes.size() << ", Archi ignorati euristica: " << edge_count << endl;
    cout << "Nuovo ordinamento salvato in ordering.txt" << endl;

    return 0;
}
