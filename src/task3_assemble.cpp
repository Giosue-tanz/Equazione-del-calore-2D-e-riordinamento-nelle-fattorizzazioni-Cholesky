#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <iomanip>
#include <sys/stat.h>

using namespace std;

// Cartella di output (stessa degli altri task)
static const string OUTPUT_DIR = "output/";

// Parametri fisici
static const double KAPPA = 0.01;

// Termine sorgente f(x,y) = exp(-10*(x^2 + y^2))
inline double f_sorgente(double x, double y) {
    return exp(-10.0 * (x * x + y * y));
}


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
    if (argc < 2 || argc > 4) {
        cerr << "\n  Assembla il sistema lineare Au = b per l'equazione del calore 2D." << endl;
        cerr << "\n  Uso:" << endl;
        cerr << "    " << argv[0] << " <N> [-r]" << endl;
        cerr << "\n  Argomenti:" << endl;
        cerr << "    N           Dimensione della griglia: genera N^2 incognite interne" << endl;
        cerr << "    -r          Applica l'ordinamento Nested Dissection (richiede ordering.txt)" << endl;
        cerr << "                Se omesso, viene usato l'ordinamento naturale (lessicografico)" << endl;
        cerr << "\n  Prerequisiti:" << endl;
        cerr << "    task1 N           deve essere stato eseguito per produrre coords.txt" << endl;
        cerr << "    task2 [--nd]      deve essere stato eseguito per produrre ordering.txt" << endl;
        cerr << "\n  Esempio:" << endl;
        cerr << "    " << argv[0] << " 64 -r" << endl << endl;
        return 1;
    }

    int N = stoi(argv[1]);
    if (N <= 0) {
        cerr << "Errore: N deve essere maggiore di 0." << endl;
        return 1;
    }

    bool usa_reorder = false;
    int scelta_bordo = 0;
    
    for (int i = 2; i < argc; ++i) {
        string flag = argv[i];
        if (flag == "-r") {
            usa_reorder = true;
        } else {
            try {
                scelta_bordo = stoi(flag);
            } catch (...) {
                cerr << "Errore: flag non riconosciuto o numero invalido '" << flag << "'." << endl;
                return 1;
            }
        }
    }

    int num_nodi = N * N;
    double h = 1.0 / (N + 1);
    double coeff_diag      = -4.0 * KAPPA / (h * h);   // entrata diagonale
    double coeff_fuori_diag = KAPPA / (h * h);          // entrata fuori diagonale

    // --------------------------------------------------------
    // Scelta della condizione al bordo
    // --------------------------------------------------------
    const string nomi_bordo[] = {
        "Omogenee          u0(x,y) = 0",
        "Costante          u0(x,y) = 10.0",
        "Gradiente lineare u0(x,y) = x",
        "Sinusoidale       u0(x,y) = sin(pi*x)",
        "Armonica          u0(x,y) = x^2 - y^2",
        "Shock termico     u0(x,y) = 1 se y=1, 0 altrove"
    };

    if (scelta_bordo < 0 || scelta_bordo > 5) {
        cerr << "\n  Errore: inserire un numero di condizione al bordo valido come argomento [0-5]." << endl;
        for (int i = 0; i <= 5; ++i) {
            cerr << "    [" << i << "]  " << nomi_bordo[i] << endl;
        }
        return 1;
    }
    
    cout << "\n  Condizioni al bordo di Dirichlet: [" << scelta_bordo << "] " << nomi_bordo[scelta_bordo] << endl;

    // --------------------------------------------------------
    // 2. Lettura di coords.txt — coordinate dei nodi interni
    //    Formato: id riga colonna x y
    // --------------------------------------------------------

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

    // --------------------------------------------------------
    // 3. Lettura di ordering.txt (solo se -r)
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

    // --------------------------------------------------------
    // 4. Assemblaggio della matrice A e del vettore b
    //    Iteriamo sui nodi nel NUOVO ordinamento (k = 0..N^2-1).
    //    Per ogni k ricaviamo l'old_id = perm[k] e da esso le
    //    coordinate geometriche (i, j, x, y).
    // --------------------------------------------------------

    // Pre-allochiamo: al più 5 entrate per nodo (diagonale + 4 vicini)
    // Il file A.txt e' simmetrico: ogni coppia (k, k') con k != k' compare due volte
    vector<Triplet> triplets;
    triplets.reserve(5 * num_nodi);

    vector<double> b(num_nodi, 0.0);

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
        b[k] = -f_sorgente(x, y);

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
                //   b[k] -= (kappa/h^2) * u0(x_bordo, y_bordo)
                double x_bordo = ni * h;
                double y_bordo = nj * h;
                b[k] -= coeff_fuori_diag * get_u0(x_bordo, y_bordo, scelta_bordo);
            } else {
                // Vicino interno: calcola il suo old_id, poi il new_id
                int old_id_vicino = id_naturale(ni, nj, N);
                int k_vicino      = inv_perm[old_id_vicino];

                // Entrata fuori diagonale (simmetrica: scritta due volte)
                triplets.push_back({k, k_vicino, coeff_fuori_diag});
            }
        }
    }

    // --------------------------------------------------------
    // 5. Scrittura di A.txt — tripli (i, j, val)
    // --------------------------------------------------------

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
    ofstream file_b(OUTPUT_DIR + "rhs.txt");
    if (!file_b.is_open()) {
        cerr << "Errore: impossibile creare " << OUTPUT_DIR << "rhs.txt" << endl;
        return 1;
    }
    file_b << fixed << setprecision(10);
    for (int k = 0; k < num_nodi; ++k) {
        file_b << b[k] << "\n";
    }
    file_b.close();

    // --------------------------------------------------------
    // 7. Riepilogo e benchmark
    // --------------------------------------------------------

    cout << "\n  Assemblaggio completato." << endl;
    cout << "  " << string(52, '-') << endl;
    cout << "    Dimensione griglia    N  = " << N << "  (" << num_nodi << " incognite)" << endl;
    cout << "    Passo spaziale        h  = " << h << endl;
    cout << "    Diffusivita' termica  k  = " << KAPPA << endl;
    cout << "    Ordinamento           :  " << (usa_reorder ? "Nested Dissection" : "Naturale (lessicografico)") << endl;
    cout << "    Condizione al bordo   :  [" << scelta_bordo << "] " << nomi_bordo[scelta_bordo] << endl;
    cout << "    Entrate scritte (A)   :  " << triplets.size() << endl;
    cout << "  " << string(52, '-') << endl;
    cout << "  File prodotti: A.txt, rhs.txt  ->  cartella output/" << endl << endl;

    return 0;
}
