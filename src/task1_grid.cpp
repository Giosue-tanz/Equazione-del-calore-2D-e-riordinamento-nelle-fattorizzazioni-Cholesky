#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>

using namespace std;

// Funzione per calcolare l'identificatore univoco del nodo a partire dagli indici (i, j)
// i, j vanno da 1 a N.
inline int getNodeId(int i, int j, int N) {
    return (i - 1) * N + (j - 1);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Uso: " << argv[0] << " <N>" << endl;
        cerr << "N: numero di nodi interni per lato (la griglia intera sara' (N+2)x(N+2))" << endl;
        return 1;
    }

    int N = stoi(argv[1]);
    if (N <= 0) {
        cerr << "Errore: N deve essere maggiore di 0." << endl;
        return 1;
    }

    int num_nodes = N * N;

    // --- Generazione Nodi (coords.txt) ---
    ofstream coords_file("coords.txt");
    if (!coords_file.is_open()) {
        cerr << "Errore: impossibile creare coords.txt" << endl;
        return 1;
    }

    // Impostiamo la precisione per la stampa delle coordinate
    coords_file << fixed << setprecision(8);

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            int n = getNodeId(i, j, N);
            double x = (double)i / (N + 1);
            double y = (double)j / (N + 1);
            
            // Formato: n i j x y
            coords_file << n << " " << i << " " << j << " " << x << " " << y << "\n";
        }
    }
    coords_file.close();

    // --- Generazione Grafo di Adiacenza (connectivity.txt) ---
    // Usiamo una lista di adiacenza come richiesto dal documento di progetto
    vector<vector<int>> adj_list(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        adj_list[i].reserve(4); // Al massimo 4 vicini per ogni nodo interno
    }

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            int u = getNodeId(i, j, N);

            // Per evitare duplicati e mantenere u < v nel file finale, 
            // basta guardare in avanti (destra e alto) per popolare gli archi non orientati.
            // Aggiungiamo comunque l'arco in entrambe le direzioni nella lista di adiacenza.
            
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

    // Scriviamo gli archi su connectivity.txt (solo u < v)
    ofstream conn_file("connectivity.txt");
    if (!conn_file.is_open()) {
        cerr << "Errore: impossibile creare connectivity.txt" << endl;
        return 1;
    }

    for (int u = 0; u < num_nodes; ++u) {
        for (int v : adj_list[u]) {
            if (u < v) {
                conn_file << u << " " << v << "\n";
            }
        }
    }
    conn_file.close();

    cout << "Griglia e grafo generati con successo per N = " << N << "." << endl;
    cout << "File creati: coords.txt, connectivity.txt" << endl;

    return 0;
}
