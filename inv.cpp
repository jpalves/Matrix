#include <iostream>
#include <vector>
#include <cmath>
#include "mat.hpp"
#include "complexo.hpp"
#include "svd.hpp"

using namespace std;

int main() {
    complexo<double> j = sqrt(-1);
    matriz<complexo<double>> A({{7, 1, 2},
                                {0, 0, 1},
                                {1, 0, 0}});

    cout << "==========================================" << endl;
    cout << "           MATRIZ ORIGINAL A              " << endl;
    cout << "==========================================" << endl;
    cout << A << endl;

    // --- 1. Decomposição LU e Inversão ---
    cout << "==========================================" << endl;
    cout << "   DECOMPOSIÇÃO LU & INVERSÃO (P*A = L*U) " << endl;
    cout << "==========================================" << endl;
    auto lu = decomposicao_lu(A);
    cout << "P:\n" << lu.P << endl;
    cout << "L:\n" << lu.L << endl;
    cout << "U:\n" << lu.U << endl;

    matriz<complexo<double>> A_inv_lu = inv_lu(A);
    cout << "Inversa de A calculada via LU:\n" << A_inv_lu << endl;
    cout << "Verificação (A * A_inv):\n" << A * A_inv_lu << endl;

    // --- 2. Decomposição SVD para Matriz Complexa (A * j) ---
    cout << "==========================================" << endl;
    cout << "   DECOMPOSIÇÃO SVD COMPLEXA (Aj = U*S*V^H)" << endl;
    cout << "==========================================" << endl;
    matriz<complexo<double>> Aj = A * j;
    cout << "Matriz Aj = A * j:\n" << Aj << endl;

    auto [U_c, S_c, Vh_c, Vt_c, Sigma_c] = svd(Aj);

    cout << "Matriz U (unitária):\n" << U_c << endl;
    cout << "Valores Singulares (S):\n";
    for (auto s : S_c) {
        cout << "  sigma = " << s << endl;
    }
    cout << "\nMatriz Sigma:\n" << Sigma_c << endl;
    cout << "Matriz V^H (unitária):\n" << Vh_c << endl;

    matriz<complexo<double>> Aj_rec = U_c * Sigma_c * Vh_c;
    cout << "Reconstrução (U * Sigma * V^H):\n" << Aj_rec << endl;

    cout << "Unitariedade U^H * U (deve ser I):\n" << U_c.hermitian() * U_c << endl;
    cout << "Unitariedade V * V^H (deve ser I):\n" << Vh_c.hermitian() * Vh_c << endl;

    cout << "Inversa de Aj via SVD (pinv):\n" << pinv(Aj) << endl;
    cout << "Verificação Aj * pinv(Aj):\n" << Aj * pinv(Aj) << endl;

    return 0;
}
