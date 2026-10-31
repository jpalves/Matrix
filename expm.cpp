// exponential matrix.
#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
#include <initializer_list>
#include "mat.h"
#include "complexo.hpp"
#include "bigfloat.hpp"

complexo<math::BigFloat> fact(int n) {
    if (n <= 1) return complexo<math::BigFloat>(1.0, 0.0);
    return n * fact(n - 1);
}

using namespace std;

int main() {
    matriz<complexo<double>> A({{0, 2, 3}, 
                                {4, 5, 6}, 
                                {7, 8, 9}});

    cout <<"expm dada pela expansão em série de Taylor"<< endl;
    cout << expm_teste(A) <<endl;
    
    cout << "exponêncial da matriz pelas aproximações de Padé" << endl;
    cout << expm(A) <<endl;
    
    return 0;
}
