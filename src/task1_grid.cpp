#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <sys/stat.h>

using namespace std;

// Cartella di output condivisa con gli altri task
static const string OUTPUT_DIR = "output/";

// Calcola l'ID naturale del nodo (i, j) con i,j in [1..N]
inline int getNodeId(int i, int j, int N) {
    return (i - 1) * N + (j - 1);
}

int main(int argc, char* argv[]) {
    int N = 0;

    if (argc >= 2) {
        try {
            N = stoi(argv[1]);
        } catch (...) {
            N = 0;
        }
    }

    while (N <= 0) {
        cout << "Inserire il numero di nodi interni per lato N (oppure digita 'q' per uscire): ";
        string input_str;
        if (!(cin >> input_str) || input_str == "q" || input_str == "Q") {
            cout << "Uscita richiesta dall'utente." << endl;
            return 0;
        }
        try {
            N = stoi(input_str);
            if (N <= 0) {
                cerr << "Errore: N deve essere un intero maggiore di 0. Riprova." << endl;
            }
        } catch (...) {
            cerr << "Errore: input non valido. Riprova." << endl;
            N = 0;
        }
    }

    int num_nodes = N * N;

    // ── Generazione Nodi (coords.txt) ────────────────────────────────────────
    ofstream coords_file(OUTPUT_DIR + "coords.txt");
    if (!coords_file.is_open()) {
        cerr << "Errore: impossibile creare " << OUTPUT_DIR << "coords.txt" << endl;
        return 1;
    }

    // Formato: id i j x y  (x = i/(N+1), y = j/(N+1))
    coords_file << fixed << setprecision(8);

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            int n  = getNodeId(i, j, N);
            double x = (double)i / (N + 1);
            double y = (double)j / (N + 1);
            coords_file << n << " " << i << " " << j << " " << x << " " << y << "\n";
        }
    }
    coords_file.close();

    // ── Generazione Grafo di Adiacenza (connectivity.txt) ────────────────────
    // Lista di adiacenza bidirezionale; nel file finale scriviamo solo u < v.
    vector<vector<int>> adj_list(num_nodes);
    for (int i = 0; i < num_nodes; ++i)
        adj_list[i].reserve(4); // al massimo 4 vicini per nodo interno

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            int u = getNodeId(i, j, N);

            // Vicino in alto (i, j+1)
            if (j < N) {
                int v_top = getNodeId(i, j + 1, N);
                adj_list[u].push_back(v_top);
                adj_list[v_top].push_back(u);
            }

            // Vicino a destra (i+1, j)
            if (i < N) {
                int v_right = getNodeId(i + 1, j, N);
                adj_list[u].push_back(v_right);
                adj_list[v_right].push_back(u);
            }
        }
    }

    ofstream conn_file(OUTPUT_DIR + "connectivity.txt");
    if (!conn_file.is_open()) {
        cerr << "Errore: impossibile creare " << OUTPUT_DIR << "connectivity.txt" << endl;
        return 1;
    }

    int edge_id = 0;
    for (int u = 0; u < num_nodes; ++u) {
        for (int v : adj_list[u]) { //Per ogni vicino v che trovo all'interno della lista dei vicini di u... 
            if (u < v) {
                conn_file << edge_id++ << " " << u << " " << v << "\n";
            }
        }
    }
    conn_file.close();

    cout << "Griglia e grafo generati con successo per N = " << N << "." << endl;
    cout << "File creati: " << OUTPUT_DIR << "coords.txt, "
         << OUTPUT_DIR << "connectivity.txt" << endl;

    return 0;
}
