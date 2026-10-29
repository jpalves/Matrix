#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "mat.hpp"

namespace mat {

// Helper para conjugação complexa
template <typename T>
constexpr T conj_scalar(const T& x) {
    if constexpr (is_complex_v<T>) {
        return std::conj(x);
    } else {
        return x;
    }
}

// Helper para construção de escalar T a partir de real e imag
template <typename T>
constexpr T make_scalar_val(double r, double i = 0.0) {
    if constexpr (is_complex_v<T>) {
        return T(r, i);
    } else {
        return static_cast<T>(r);
    }
}

// Estrutura que armazena a decomposição SVD: A = U * Sigma * V^H
template <Scalar T>
struct SVDResult {
    Matrix<T> U;                        // Matriz unitária/ortogonal m x m
    std::vector<real_type_t<T>> S;      // Valores singulares ordenados (sigma_0 >= sigma_1 >= ...)
    Matrix<T> Vh;                       // Matriz V conjugada transposta (V^H) n x n
    Matrix<T> Vt;                       // Alias para Vh (mantém retrocompatibilidade)
    Matrix<T> Sigma;                    // Matriz diagonal m x n contendo os valores singulares
};

// Algoritmo SVD de Jacobi Unilateral Geral (Real e Complexo / Unitary Hestenes)
// Suporta números reais (float, double) e complexos (std::complex, complexo<T>).
// Suporta matrizes quadradas e retangulares (m >= n e m < n).
template <Scalar T>
SVDResult<T> svd(const Matrix<T>& A, double tol = 1e-15, std::size_t max_sweeps = 100) {
    using real_t = real_type_t<T>;
    const std::size_t m = A.rows();
    const std::size_t n = A.cols();

    if (m == 0 || n == 0) {
        throw std::invalid_argument("A matriz não pode ser vazia para a decomposição SVD.");
    }

    // Se m < n, decompomos A^H: A^H = U_t * Sigma_t * V_t^H => A = V_t * Sigma_t^T * U_t^H
    if (m < n) {
        auto res_t = svd(A.hermitian(), tol, max_sweeps);
        Matrix<T> Sigma(m, n, T{0});
        for (std::size_t i = 0; i < res_t.S.size(); ++i) {
            Sigma(i, i) = make_scalar_val<T>(res_t.S[i]);
        }
        Matrix<T> U = res_t.Vh.hermitian();
        Matrix<T> Vh = res_t.U.hermitian();
        return {U, res_t.S, Vh, Vh, Sigma};
    }

    // Caso m >= n
    Matrix<T> B = A; // B convergirá para U * Sigma (colunas mutuamente ortogonais)
    Matrix<T> V = Matrix<T>::identity(n);

    for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep) {
        double max_off_diag = 0.0;

        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = j + 1; k < n; ++k) {
                double alpha = 0.0;
                double beta  = 0.0;
                std::complex<double> gamma(0.0, 0.0);

                for (std::size_t i = 0; i < m; ++i) {
                    std::complex<double> bj = B(i, j);
                    std::complex<double> bk = B(i, k);
                    alpha += std::norm(bj);
                    beta  += std::norm(bk);
                    gamma += std::conj(bj) * bk;
                }

                double gamma_abs = std::abs(gamma);
                double denom = std::sqrt(alpha * beta) + 1e-30;
                double err = gamma_abs / denom;
                if (err > max_off_diag) {
                    max_off_diag = err;
                }

                if (gamma_abs <= tol * std::sqrt(alpha * beta)) {
                    continue;
                }

                // Rotação de Jacobi complexa (Unitária 2x2)
                double zeta = (alpha - beta) / (2.0 * gamma_abs);
                double t = (zeta >= 0.0) ? (1.0 / (zeta + std::sqrt(1.0 + zeta * zeta)))
                                         : (-1.0 / (-zeta + std::sqrt(1.0 + zeta * zeta)));
                double c = 1.0 / std::sqrt(1.0 + t * t);
                double s0 = t * c;

                std::complex<double> s_c = s0 * (std::conj(gamma) / gamma_abs);
                std::complex<double> s_bar = std::conj(s_c);

                T c_val = make_scalar_val<T>(c);
                T s_val = make_scalar_val<T>(s_c.real(), s_c.imag());
                T s_bar_val = make_scalar_val<T>(s_bar.real(), s_bar.imag());

                // Aplica rotação unitária às colunas j e k de B
                for (std::size_t i = 0; i < m; ++i) {
                    T b_ij = B(i, j);
                    T b_ik = B(i, k);
                    B(i, j) = c_val * b_ij + s_val * b_ik;
                    B(i, k) = -s_bar_val * b_ij + c_val * b_ik;
                }

                // Acumula rotação unitária em V
                for (std::size_t i = 0; i < n; ++i) {
                    T v_ij = V(i, j);
                    T v_ik = V(i, k);
                    V(i, j) = c_val * v_ij + s_val * v_ik;
                    V(i, k) = -s_bar_val * v_ij + c_val * v_ik;
                }
            }
        }

        if (max_off_diag < tol) {
            break;
        }
    }

    // Extrai valores singulares (normas euclidianas das colunas de B) e normaliza U_thin
    std::vector<real_t> S(n, 0.0);
    Matrix<T> U_thin(m, n, T{0});

    for (std::size_t j = 0; j < n; ++j) {
        double norm = 0.0;
        for (std::size_t i = 0; i < m; ++i) {
            norm += std::norm(static_cast<std::complex<double>>(B(i, j)));
        }
        norm = std::sqrt(norm);
        S[j] = static_cast<real_t>(norm);

        if (norm > 1e-15) {
            for (std::size_t i = 0; i < m; ++i) {
                U_thin(i, j) = B(i, j) / make_scalar_val<T>(norm);
            }
        }
    }

    // Ordenação decrescente dos valores singulares
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&S](std::size_t a, std::size_t b) {
        return S[a] > S[b];
    });

    std::vector<real_t> S_sorted(n);
    Matrix<T> U_sorted(m, n, T{0});
    Matrix<T> V_sorted(n, n, T{0});

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t src = order[col];
        S_sorted[col] = S[src];
        for (std::size_t i = 0; i < m; ++i) {
            U_sorted(i, col) = U_thin(i, src);
        }
        for (std::size_t i = 0; i < n; ++i) {
            V_sorted(i, col) = V(i, src);
        }
    }

    // Matriz Sigma (m x n)
    Matrix<T> Sigma(m, n, T{0});
    for (std::size_t i = 0; i < n; ++i) {
        Sigma(i, i) = make_scalar_val<T>(S_sorted[i]);
    }

    // Completa base ortonormal para U (m x m) caso m > n (Full SVD)
    Matrix<T> U_full(m, m, T{0});
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            U_full(i, j) = U_sorted(i, j);
        }
    }
    if (m > n) {
        for (std::size_t j = n; j < m; ++j) {
            std::vector<T> col(m, T{0});
            col[j] = T{1};
            for (std::size_t k = 0; k < j; ++k) {
                T dot = T{0};
                for (std::size_t i = 0; i < m; ++i) {
                    dot += conj_scalar(U_full(i, k)) * col[i];
                }
                for (std::size_t i = 0; i < m; ++i) {
                    col[i] -= dot * U_full(i, k);
                }
            }
            double norm = 0.0;
            for (std::size_t i = 0; i < m; ++i) {
                norm += std::norm(static_cast<std::complex<double>>(col[i]));
            }
            norm = std::sqrt(norm);
            if (norm > 1e-12) {
                for (std::size_t i = 0; i < m; ++i) {
                    U_full(i, j) = col[i] / make_scalar_val<T>(norm);
                }
            }
        }
    }

    Matrix<T> Vh = V_sorted.hermitian().chop(1e-14);
    U_full.chop(1e-14);

    return {U_full, S_sorted, Vh, Vh, Sigma};
}

// Cria uma matriz diagonal a partir de um std::vector
template <Scalar T, typename V>
Matrix<T> diag(const std::vector<V>& v) {
    std::size_t n = v.size();
    Matrix<T> res(n, n, T{0});
    for (std::size_t i = 0; i < n; ++i) {
        res(i, i) = make_scalar_val<T>(v[i]);
    }
    return res;
}

// Inversão / Pseudo-inversa via SVD: A^+ = V * Sigma^+ * U^H
template <Scalar T>
Matrix<T> pinv(const Matrix<T>& A, double rcond = 1e-15, bool refine = true) {
    auto [U, S, Vh, Vt, Sigma] = svd(A);
    const std::size_t m = A.rows();
    const std::size_t n = A.cols();

    Matrix<T> Sigma_plus(n, m, T{0});
    real_type_t<T> max_s = S.empty() ? 0 : S[0];
    real_type_t<T> cutoff = rcond * max_s;

    for (std::size_t i = 0; i < S.size(); ++i) {
        if (S[i] > cutoff) {
            Sigma_plus(i, i) = make_scalar_val<T>(1.0 / S[i]);
        }
    }

    // A^+ = V * Sigma_plus * U^H (como Vh = V^H, temos V = Vh^H)
    Matrix<T> X = Vh.hermitian() * Sigma_plus * U.hermitian();

    // Refinamento iterativo de Newton-Schultz: eleva a precisão ao limite da máquina
    if (refine && A.is_square()) {
        Matrix<T> I = Matrix<T>::identity(n);
        Matrix<T> E = I - A * X;
        X += X * E;
    }

    return X.chop(1e-14);
}

template <Scalar T>
Matrix<T> inv_svd(const Matrix<T>& A, double rcond = 1e-15, bool refine = true) {
    return pinv(A, rcond, refine);
}

} // namespace mat

// Exporta para escopo global para compatibilidade direta
using mat::svd;
using mat::SVDResult;
using mat::diag;
using mat::pinv;
using mat::inv_svd;
