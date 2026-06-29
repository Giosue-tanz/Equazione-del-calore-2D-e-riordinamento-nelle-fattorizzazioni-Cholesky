// ============================================================
//  task3_assemble.cpp
//  Modulo 3 — Assemblaggio del sistema lineare Au = b
//
//  Equazione del calore 2D stazionaria con differenze finite:
//      kappa * (u_xx + u_yy) + f(x,y) = 0,  (x,y) in Omega=[0,1]^2
//      u = 0 sul bordo
//
//  Parametri fisici:
//      kappa = 0.01
//      f(x,y) = exp(-10*(x^2 + y^2))
//
//  Stencil a 5 punti:
//      A[k,k]   = -4*kappa/h^2   (diagonale)
//      A[k,k']  = +kappa/h^2     (fuori diagonale, per ogni vicino interno)
//      b[k]     = -f(x_i, y_j)   (con eventuale correzione Dirichlet, qui nulla)
//
//  Uso:
//      ./task3_assemble <N>
//      ./task3_assemble <N> --reorder
//
//  Output:
//      output/A.txt   — tripli (i, j, val) per ogni entrata non nulla (simmetrica)
//      output/rhs.txt — vettore dei termini noti, uno per riga
// ============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <chrono>
#include <sys/stat.h>

using namespace std;
using namespace chrono;

// Cartella di output (stessa degli altri task)
static const string OUTPUT_DIR = "output/";

// Parametri fisici
static const double KAPPA = 0.01;

// Termine sorgente f(x,y) = exp(-10*(x^2 + y^2))
inline double f_sorgente(double x, double y) {
    return exp(-10.0 * (x * x + y * y));
}

// ============================================================
//  Parte Opzionale — Condizioni al bordo non omogenee u0(x,y)
//
//  Scelta 0: u0 = 0              (omogenee, default)
//  Scelta 1: u0 = 10.0           (costante — sanity check)
//  Scelta 2: u0 = x              (gradiente lineare)
//  Scelta 3: u0 = sin(pi*x)      (sinusoidale continua)
//  Scelta 4: u0 = x^2 - y^2     (armonica, Laplaciano = 0)
//  Scelta 5: u0 = 1 se y=1, 0   (shock termico bordo superiore)
// ============================================================
double get_u0(double x, double y, int scelta) {
    switch (scelta) {
        case 0: return 0.0;
        case 1: return 10.0;
        case 2: return x;
        case 3: return sin(M_PI * x);
        case 4: return x * x - y * y;
        case 5: return (y >= 1.0 - 1e-10) ? 1.0 : 0.0;
        default:
            cerr << "Errore: scelta bordo non valida (" << scelta << ")." << endl;
            exit(1);
    }
}

// Calcola l'ID naturale del nodo (i, j) con i,j in [1..N]
// ID = (i-1)*N + (j-1)  => riga maggiore, colonna minore
inline int id_naturale(int i, int j, int N) {
    return (i - 1) * N + (j - 1);
}

// Struttura per i tripli (i, j, val) della matrice sparsa
struct Triplet {
    int row;
    int col;
    double val;
};

int main(int argc, char* argv[]) {
    // --------------------------------------------------------
    // 1. Parsing degli argomenti da riga di comando
    // --------------------------------------------------------
    if (argc < 2 || argc > 3) {
        cerr << "Uso: " << argv[0] << " <N> [--reorder]" << endl;
        cerr << "  N         : numero di nodi interni per lato" << endl;
        cerr << "  --reorder : usa l'ordinamento nested dissection da ordering.txt" << endl;
        return 1;
    }

    int N = stoi(argv[1]);
    if (N <= 0) {
        cerr << "Errore: N deve essere maggiore di 0." << endl;
        return 1;
    }

    bool usa_reorder = false;
    if (argc == 3) {
        string flag = argv[2];
        if (flag == "--reorder") {
            usa_reorder = true;
        } else {
            cerr << "Errore: flag non riconosciuto '" << flag << "'. Usa --reorder." << endl;
            return 1;
        }
    }

    int num_nodi = N * N;
    double h = 1.0 / (N + 1);
    double coeff_diag      = -4.0 * KAPPA / (h * h);   // entrata diagonale
    double coeff_fuori_diag = KAPPA / (h * h);          // entrata fuori diagonale

    // --------------------------------------------------------
    // Menu interattivo: scelta della condizione al bordo
    // --------------------------------------------------------
    const string nomi_bordo[] = {
        "Omogenee          u0(x,y) = 0",
        "Costante          u0(x,y) = 10.0",
        "Gradiente lineare u0(x,y) = x",
        "Sinusoidale       u0(x,y) = sin(pi*x)",
        "Armonica          u0(x,y) = x^2 - y^2",
        "Shock termico     u0(x,y) = 1 se y=1, 0 altrove"
    };
    int scelta_bordo = 0;
    cout << "\nCondizioni al bordo di Dirichlet:" << endl;
    for (int i = 0; i <= 5; ++i)
        cout << "  [" << i << "] " << nomi_bordo[i] << endl;
    cout << "Scelta [0-5]: ";
    if (!(cin >> scelta_bordo) || scelta_bordo < 0 || scelta_bordo > 5) {
        cerr << "Errore: inserire un valore intero tra 0 e 5." << endl;
        return 1;
    }

    mkdir(OUTPUT_DIR.c_str(), 0755);

    // --------------------------------------------------------
    // 2. Lettura di coords.txt — coordinate dei nodi interni
    //    Formato: id riga colonna x y
    // --------------------------------------------------------
    auto t_inizio_lettura = high_resolution_clock::now();

    struct NodeCoord {
        int id;
        int i;
        int j;
        double x;
        double y;
    };

    vector<NodeCoord> coords(num_nodi);

    ifstream file_coords(OUTPUT_DIR + "coords.txt");
    if (!file_coords.is_open()) {
        cerr << "Errore: impossibile aprire " << OUTPUT_DIR << "coords.txt" << endl;
        cerr << "Assicurarsi di aver eseguito task1_grid prima." << endl;
        return 1;
    }
    {
        int id, ri, rj;
        double rx, ry;
        int contatore_nodi = 0;
        while (file_coords >> id >> ri >> rj >> rx >> ry) {
            if (id < 0 || id >= num_nodi) {
                cerr << "Errore: ID nodo " << id << " fuori range [0, " << num_nodi - 1 << "]" << endl;
                return 1;
            }
            coords[id] = {id, ri, rj, rx, ry};
            contatore_nodi++;
        }
        if (contatore_nodi != num_nodi) {
            cerr << "Errore: lette " << contatore_nodi << " righe da coords.txt, ma N=" << N << " richiede " << num_nodi << " nodi." << endl;
            cerr << "Assicurati di aver generato coords.txt con lo stesso N." << endl;
            return 1;
        }
    }
    file_coords.close();

    auto t_fine_lettura_coords = high_resolution_clock::now();

    // --------------------------------------------------------
    // 3. Lettura di ordering.txt (solo se --reorder)
    //    Formato: new_id old_id
    //    inv_perm[old_id] = new_id  (da ordinamento naturale a riordinato)
    // --------------------------------------------------------

    // perm[new_id]     = old_id  (da posizione riordinata a nodo originale)
    // inv_perm[old_id] = new_id  (da nodo originale a posizione riordinata)
    vector<int> perm(num_nodi);
    vector<int> inv_perm(num_nodi);

    if (usa_reorder) {
        ifstream file_ordering(OUTPUT_DIR + "ordering.txt");
        if (!file_ordering.is_open()) {
            cerr << "Errore: impossibile aprire " << OUTPUT_DIR << "ordering.txt" << endl;
            cerr << "Assicurarsi di aver eseguito task2_reorder prima." << endl;
            return 1;
        }
        int new_id, old_id;
        int contatore_ord = 0;
        while (file_ordering >> new_id >> old_id) {
            if (new_id < 0 || new_id >= num_nodi || old_id < 0 || old_id >= num_nodi) {
                cerr << "Errore: indice fuori range in ordering.txt" << endl;
                return 1;
            }
            perm[new_id]     = old_id;
            inv_perm[old_id] = new_id;
            contatore_ord++;
        }
        if (contatore_ord != num_nodi) {
            cerr << "Errore: lette " << contatore_ord << " righe da ordering.txt, ma N=" << N << " richiede " << num_nodi << " entry." << endl;
            cerr << "Assicurati di aver generato ordering.txt con lo stesso N." << endl;
            return 1;
        }
        file_ordering.close();
    } else {
        // Ordinamento naturale: identità
        for (int k = 0; k < num_nodi; ++k) {
            perm[k]     = k;
            inv_perm[k] = k;
        }
    }

    auto t_fine_lettura_ordering = high_resolution_clock::now();

    // --------------------------------------------------------
    // 4. Assemblaggio della matrice A e del vettore b
    //    Iteriamo sui nodi nel NUOVO ordinamento (k = 0..N^2-1).
    //    Per ogni k ricaviamo l'old_id = perm[k] e da esso le
    //    coordinate geometriche (i, j, x, y).
    // --------------------------------------------------------
    auto t_inizio_assemblaggio = high_resolution_clock::now();

    // Pre-allochiamo: al più 5 entrate per nodo (diagonale + 4 vicini)
    // Il file A.txt e' simmetrico: ogni coppia (k, k') con k != k' compare due volte
    vector<Triplet> triplets;
    triplets.reserve(5 * num_nodi);

    vector<double> rhs(num_nodi, 0.0);

    for (int k = 0; k < num_nodi; ++k) {
        int old_id = perm[k];           // nodo nel sistema naturale
        const NodeCoord& nc = coords[old_id];
        int ri = nc.i;
        int rj = nc.j;
        double x = nc.x;
        double y = nc.y;

        // --- Entrata diagonale ---
        triplets.push_back({k, k, coeff_diag});

        // --- Termine noto: b_k = -f(x_i, y_j) ---
        // (il segno meno viene dal fatto che lo stencil e' negativo e noi assembliamo Au=b)
        rhs[k] = -f_sorgente(x, y);

        // --- Vicini e fuori diagonale ---
        // Direzioni: sinistra (i-1,j), destra (i+1,j), basso (i,j-1), alto (i,j+1)
        const int di[] = {-1, +1,  0,  0};
        const int dj[] = { 0,  0, -1, +1};

        for (int d = 0; d < 4; ++d) {
            int ni = ri + di[d];
            int nj = rj + dj[d];

            if (ni < 1 || ni > N || nj < 1 || nj > N) {
                // Il vicino cade sul bordo: il suo valore e' noto (condizione di Dirichlet).
                // Non entra nella matrice A, ma contribuisce al termine noto:
                //   rhs[k] -= (kappa/h^2) * u0(x_bordo, y_bordo)
                double x_bordo = ni * h;
                double y_bordo = nj * h;
                rhs[k] -= coeff_fuori_diag * get_u0(x_bordo, y_bordo, scelta_bordo);
            } else {
                // Vicino interno: calcola il suo old_id, poi il new_id
                int old_id_vicino = id_naturale(ni, nj, N);
                int k_vicino      = inv_perm[old_id_vicino];

                // Entrata fuori diagonale (simmetrica: scritta due volte)
                triplets.push_back({k, k_vicino, coeff_fuori_diag});
            }
        }
    }

    auto t_fine_assemblaggio = high_resolution_clock::now();

    // --------------------------------------------------------
    // 5. Scrittura di A.txt — tripli (i, j, val)
    // --------------------------------------------------------
    auto t_inizio_scrittura = high_resolution_clock::now();

    ofstream file_A(OUTPUT_DIR + "A.txt");
    if (!file_A.is_open()) {
        cerr << "Errore: impossibile creare " << OUTPUT_DIR << "A.txt" << endl;
        return 1;
    }
    file_A << fixed << setprecision(10);
    for (const auto& t : triplets) {
        file_A << t.row << " " << t.col << " " << t.val << "\n";
    }
    file_A.close();

    // --------------------------------------------------------
    // 6. Scrittura di rhs.txt — un valore per riga
    // --------------------------------------------------------
    ofstream file_rhs(OUTPUT_DIR + "rhs.txt");
    if (!file_rhs.is_open()) {
        cerr << "Errore: impossibile creare " << OUTPUT_DIR << "rhs.txt" << endl;
        return 1;
    }
    file_rhs << fixed << setprecision(10);
    for (int k = 0; k < num_nodi; ++k) {
        file_rhs << rhs[k] << "\n";
    }
    file_rhs.close();

    auto t_fine_scrittura = high_resolution_clock::now();

    // --------------------------------------------------------
    // 7. Riepilogo e benchmark
    // --------------------------------------------------------
    duration<double, milli> tempo_coords    = t_fine_lettura_coords    - t_inizio_lettura;
    duration<double, milli> tempo_ordering  = t_fine_lettura_ordering  - t_fine_lettura_coords;
    duration<double, milli> tempo_assembl   = t_fine_assemblaggio      - t_inizio_assemblaggio;
    duration<double, milli> tempo_scrittura = t_fine_scrittura         - t_inizio_scrittura;
    duration<double, milli> tempo_totale    = t_fine_scrittura         - t_inizio_lettura;

    cout << "Assemblaggio completato con successo." << endl;
    cout << "  N                  : " << N << endl;
    cout << "  Nodi interni (N^2) : " << num_nodi << endl;
    cout << "  Passo h            : " << h << endl;
    cout << "  kappa              : " << KAPPA << endl;
    cout << "  Ordinamento        : " << (usa_reorder ? "Nested Dissection (--reorder)" : "Naturale") << endl;
    cout << "  Condizione bordo   : [" << scelta_bordo << "] " << nomi_bordo[scelta_bordo] << endl;
    cout << "  Entrate in A.txt   : " << triplets.size()
         << "  (attese: " << num_nodi + 4LL * num_nodi << " max)" << endl;
    cout << "\n=== TEMPI DI ESECUZIONE (BENCHMARK) ===" << endl;
    cout << fixed << setprecision(3);
    cout << "  Lettura coords.txt   : " << tempo_coords.count()    << " ms" << endl;
    if (usa_reorder)
    cout << "  Lettura ordering.txt : " << tempo_ordering.count()  << " ms" << endl;
    cout << "  Assemblaggio matrice : " << tempo_assembl.count()   << " ms" << endl;
    cout << "  Scrittura A+rhs.txt  : " << tempo_scrittura.count() << " ms" << endl;
    cout << "  Tempo Totale         : " << tempo_totale.count()    << " ms" << endl;

    return 0;
}
