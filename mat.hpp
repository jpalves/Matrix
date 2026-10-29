// mat.hpp (e mat.h)
//
// Biblioteca de matrizes densas em C++20 STL moderno.
// Reimplementação moderna baseada no estilo de shipley_coleman.hpp:
// - Gerenciamento automático de memória (std::vector, Rule of Zero).
// - Compatibilidade com std::span, std::ranges e algoritmos STL (<algorithm>, <numeric>).
// - Suporte a tipos reais (floating_point) e complexos (std::complex).
// - Indexação bidimensional funcional (i, j), subscript m[i][j] via std::span e C++23 m[i, j].
// - Operações de álgebra linear: determinante (LU com pivotação parcial),
//   inversão nodal in situ (Shipley-Coleman com pivotação parcial),
//   inversão por adjunta/cofatores (Cramer), transposição e operadores aritméticos.
// - 100% de retrocompatibilidade com a API clássica do mat.h (matriz<T>, linhas, colunas, inv, Inv, etc.).

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <ostream>
#include <ranges>
#include <span>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace mat {

// Verificação de tipos complexos (std::complex ou classes derivadas como complexo<T>)
template <typename T>
struct is_complex : std::bool_constant<
    requires {
      typename T::value_type;
      requires (std::same_as<T, std::complex<typename T::value_type>> ||
                std::derived_from<T, std::complex<typename T::value_type>>);
    }> {};

template <typename T>
inline constexpr bool is_complex_v = is_complex<T>::value;

// Traits para extração do tipo escalar real
template <typename T, typename = void>
struct scalar_traits {
  using real_type = T;
};

template <typename T>
struct scalar_traits<T, std::enable_if_t<is_complex_v<T>>> {
  using real_type = typename T::value_type;
};

template <typename T>
using real_type_t = typename scalar_traits<T>::real_type;

// Conceito moderno C++20 para tipos escalares suportados
template <typename T>
concept Scalar = std::floating_point<T> ||
                 (is_complex_v<T> && std::floating_point<typename T::value_type>);

// Classe Matrix: matriz densa em linha contígua (row-major)
template <Scalar T>
class Matrix {
 public:
  using value_type      = T;
  using size_type       = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference       = T&;
  using const_reference = const T&;
  using pointer         = T*;
  using const_pointer   = const T*;
  using iterator        = typename std::vector<T>::iterator;
  using const_iterator  = typename std::vector<T>::const_iterator;

  // --- Construtores ---
  Matrix() = default;

  Matrix(std::size_t rows, std::size_t cols, const T& init = T{})
      : rows_(rows), cols_(cols), data_(rows * cols, init) {}

  explicit Matrix(std::size_t n, const T& init = T{}) : Matrix(n, n, init) {}

  Matrix(std::size_t rows, std::size_t cols, const T* arr)
      : rows_(rows), cols_(cols), data_(arr, arr + (rows * cols)) {}

  Matrix(std::size_t rows, std::size_t cols, std::span<const T> sp)
      : rows_(rows), cols_(cols), data_(sp.begin(), sp.end()) {
    if (sp.size() != rows * cols) {
      throw std::invalid_argument("Tamanho do span incompatível com as dimensões da matriz");
    }
  }

  Matrix(std::initializer_list<std::initializer_list<T>> list) {
    rows_ = list.size();
    cols_ = (rows_ > 0) ? list.begin()->size() : 0;
    data_.reserve(rows_ * cols_);
    for (const auto& row_list : list) {
      if (row_list.size() != cols_) {
        throw std::invalid_argument("Todas as linhas devem ter a mesma dimensão");
      }
      data_.insert(data_.end(), row_list.begin(), row_list.end());
    }
  }

  // --- Acesso a elementos ---
  [[nodiscard]] constexpr T& operator()(std::size_t i, std::size_t j) noexcept {
    return data_[i * cols_ + j];
  }

  [[nodiscard]] constexpr const T& operator()(std::size_t i, std::size_t j) const noexcept {
    return data_[i * cols_ + j];
  }

  // Acesso via operador [] retornando span de linha (permite sintaxe m[i][j])
  [[nodiscard]] constexpr std::span<T> operator[](std::size_t i) noexcept {
    return row(i);
  }

  [[nodiscard]] constexpr std::span<const T> operator[](std::size_t i) const noexcept {
    return row(i);
  }

#if defined(__cpp_multidimensional_subscript) && __cpp_multidimensional_subscript >= 202110L
  [[nodiscard]] constexpr T& operator[](std::size_t i, std::size_t j) noexcept {
    return (*this)(i, j);
  }

  [[nodiscard]] constexpr const T& operator[](std::size_t i, std::size_t j) const noexcept {
    return (*this)(i, j);
  }
#endif

  // Visão em span de cada linha contígua
  [[nodiscard]] constexpr std::span<T> row(std::size_t i) noexcept {
    return std::span<T>(data_.data() + i * cols_, cols_);
  }

  [[nodiscard]] constexpr std::span<const T> row(std::size_t i) const noexcept {
    return std::span<const T>(data_.data() + i * cols_, cols_);
  }

  // Dimensões e estado
  [[nodiscard]] constexpr std::size_t rows() const noexcept { return rows_; }
  [[nodiscard]] constexpr std::size_t cols() const noexcept { return cols_; }
  [[nodiscard]] constexpr std::size_t linhas() const noexcept { return rows_; }
  [[nodiscard]] constexpr std::size_t colunas() const noexcept { return cols_; }
  [[nodiscard]] constexpr std::size_t size() const noexcept { return data_.size(); }
  [[nodiscard]] constexpr std::size_t numel() const noexcept { return data_.size(); }
  [[nodiscard]] constexpr bool empty() const noexcept { return data_.empty(); }
  [[nodiscard]] constexpr bool is_square() const noexcept { return rows_ == cols_; }

  // Acesso contíguo
  [[nodiscard]] constexpr T* data() noexcept { return data_.data(); }
  [[nodiscard]] constexpr const T* data() const noexcept { return data_.data(); }
  [[nodiscard]] constexpr T* as_vector() noexcept { return data_.data(); }
  [[nodiscard]] constexpr const T* as_vector() const noexcept { return data_.data(); }
  [[nodiscard]] constexpr T* asVector() noexcept { return data_.data(); }
  [[nodiscard]] constexpr const T* asVector() const noexcept { return data_.data(); }

  [[nodiscard]] constexpr std::span<T> as_span() noexcept { return std::span<T>(data_); }
  [[nodiscard]] constexpr std::span<const T> as_span() const noexcept {
    return std::span<const T>(data_);
  }

  // Iteradores STL (para std::ranges e algoritmos STL)
  [[nodiscard]] constexpr iterator begin() noexcept { return data_.begin(); }
  [[nodiscard]] constexpr iterator end() noexcept { return data_.end(); }
  [[nodiscard]] constexpr const_iterator begin() const noexcept { return data_.begin(); }
  [[nodiscard]] constexpr const_iterator end() const noexcept { return data_.end(); }
  [[nodiscard]] constexpr const_iterator cbegin() const noexcept { return data_.cbegin(); }
  [[nodiscard]] constexpr const_iterator cend() const noexcept { return data_.cend(); }

  // --- Operadores Compostos ---
  Matrix<T>& operator+=(const Matrix<T>& rhs) {
    if (rows_ != rhs.rows_ || cols_ != rhs.cols_) {
      throw std::invalid_argument("Dimensões incompatíveis para adição");
    }
    std::transform(data_.begin(), data_.end(), rhs.data_.begin(), data_.begin(), std::plus<>{});
    return *this;
  }

  Matrix<T>& operator-=(const Matrix<T>& rhs) {
    if (rows_ != rhs.rows_ || cols_ != rhs.cols_) {
      throw std::invalid_argument("Dimensões incompatíveis para subtração");
    }
    std::transform(data_.begin(), data_.end(), rhs.data_.begin(), data_.begin(), std::minus<>{});
    return *this;
  }

  Matrix<T>& operator+=(const T& s) {
    for (auto& x : data_) {
      x += s;
    }
    return *this;
  }

  Matrix<T>& operator-=(const T& s) {
    for (auto& x : data_) {
      x -= s;
    }
    return *this;
  }

  Matrix<T>& operator*=(const T& s) {
    std::ranges::transform(data_, data_.begin(), [s](const T& x) { return x * s; });
    return *this;
  }

  Matrix<T>& operator/=(const T& s) {
    std::ranges::transform(data_, data_.begin(), [s](const T& x) { return x / s; });
    return *this;
  }

  // --- Permutações de Linhas e Colunas ---
  void swap_rows(std::size_t i, std::size_t j) {
    if (i >= rows_ || j >= rows_) {
      throw std::out_of_range("Índice de linha fora dos limites");
    }
    if (i == j) return;
    auto ri = row(i);
    auto rj = row(j);
    std::swap_ranges(ri.begin(), ri.end(), rj.begin());
  }

  void trocaLinhas(std::size_t i, std::size_t j) { swap_rows(i, j); }

  void swap_cols(std::size_t i, std::size_t j) {
    if (i >= cols_ || j >= cols_) {
      throw std::out_of_range("Índice de coluna fora dos limites");
    }
    if (i == j) return;
    for (std::size_t r = 0; r < rows_; ++r) {
      std::swap((*this)(r, i), (*this)(r, j));
    }
  }

  void trocaColunas(std::size_t i, std::size_t j) { swap_cols(i, j); }

  // --- Transposição e Identidade ---
  [[nodiscard]] Matrix<T> transpose() const {
    Matrix<T> res(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
      for (std::size_t j = 0; j < cols_; ++j) {
        res(j, i) = (*this)(i, j);
      }
    }
    return res;
  }

  [[nodiscard]] Matrix<T> transposta() const { return transpose(); }

  // --- Conjugação e Transposição Conjugada (Hermitiana) ---
  [[nodiscard]] Matrix<T> conjugate() const {
    Matrix<T> res = *this;
    for (auto& x : res.data_) {
      if constexpr (is_complex_v<T>) {
        x = std::conj(x);
      }
    }
    return res;
  }

  [[nodiscard]] Matrix<T> conjugada() const { return conjugate(); }

  [[nodiscard]] Matrix<T> hermitian() const {
    Matrix<T> res(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
      for (std::size_t j = 0; j < cols_; ++j) {
        if constexpr (is_complex_v<T>) {
          res(j, i) = std::conj((*this)(i, j));
        } else {
          res(j, i) = (*this)(i, j);
        }
      }
    }
    return res;
  }

  [[nodiscard]] Matrix<T> adjoint() const { return hermitian(); }
  [[nodiscard]] Matrix<T> H() const { return hermitian(); }
  [[nodiscard]] Matrix<T> transposta_conjugada() const { return hermitian(); }

  [[nodiscard]] static Matrix<T> identity(std::size_t n) {
    Matrix<T> res(n, n, T{0});
    for (std::size_t i = 0; i < n; ++i) {
      res(i, i) = T{1};
    }
    return res;
  }

  [[nodiscard]] static Matrix<T> eye(std::size_t n) { return identity(n); }

  [[nodiscard]] Matrix<T> eye() const {
    if (!is_square()) {
      throw std::invalid_argument("eye() sem argumentos requer matriz quadrada");
    }
    return identity(rows_);
  }

  // --- Limpeza de ruído numérico de ponto flutuante ---
  Matrix<T>& chop(real_type_t<T> eps = 1e-14) {
    for (auto& x : data_) {
      if (std::abs(x) < eps) {
        x = T{0};
      } else if constexpr (is_complex_v<T>) {
        using r_t = real_type_t<T>;
        r_t r = std::real(x);
        r_t im = std::imag(x);
        if (std::abs(r) < eps) r = 0;
        if (std::abs(im) < eps) im = 0;
        x = T(r, im);
      }
    }
    return *this;
  }

  [[nodiscard]] Matrix<T> chopped(real_type_t<T> eps = 1e-14) const {
    Matrix<T> res = *this;
    return res.chop(eps);
  }

  Matrix<T>& limpa(real_type_t<T> eps = 1e-14) { return chop(eps); }
  [[nodiscard]] Matrix<T> limpa(real_type_t<T> eps = 1e-14) const { return chopped(eps); }


  // --- Submatriz / Menor Complementar (Cofator) ---
  [[nodiscard]] Matrix<T> submatrix(std::size_t row_to_remove, std::size_t col_to_remove) const {
    if (rows_ <= 1 || cols_ <= 1) {
      throw std::logic_error("Matriz demasiado pequena para extrair submatriz");
    }
    if (row_to_remove >= rows_ || col_to_remove >= cols_) {
      throw std::out_of_range("Índice de linha ou coluna fora dos limites em submatrix");
    }
    Matrix<T> res(rows_ - 1, cols_ - 1);
    std::size_t r_dst = 0;
    for (std::size_t r = 0; r < rows_; ++r) {
      if (r == row_to_remove) continue;
      std::size_t c_dst = 0;
      for (std::size_t c = 0; c < cols_; ++c) {
        if (c == col_to_remove) continue;
        res(r_dst, c_dst) = (*this)(r, c);
        ++c_dst;
      }
      ++r_dst;
    }
    return res;
  }

  [[nodiscard]] Matrix<T> coff(std::size_t i, std::size_t j) const {
    return submatrix(i, j);
  }

  // --- Determinante (Eliminação Gaussiana com Pivotação Parcial) ---
  [[nodiscard]] T det() const {
    if (!is_square()) {
      throw std::invalid_argument("A matriz não é quadrada para cálculo de determinante");
    }
    const std::size_t n = rows_;
    if (n == 0) return T{1};
    if (n == 1) return (*this)(0, 0);
    if (n == 2) {
      return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
    }

    Matrix<T> temp = *this;
    T det_val = T{1};
    using real_t = real_type_t<T>;
    using std::abs;
    const real_t eps = std::numeric_limits<real_t>::epsilon() * 100;

    for (std::size_t k = 0; k < n; ++k) {
      // Procura pivô máximo na coluna k (de k até n-1)
      std::size_t pivot_row = k;
      real_t max_val = abs(temp(k, k));
      for (std::size_t i = k + 1; i < n; ++i) {
        real_t val = abs(temp(i, k));
        if (val > max_val) {
          max_val = val;
          pivot_row = i;
        }
      }

      if (max_val < eps) {
        return T{0};  // Singular
      }

      if (pivot_row != k) {
        temp.swap_rows(k, pivot_row);
        det_val = -det_val;
      }

      const T pivot = temp(k, k);
      det_val *= pivot;

      const auto rk = temp.row(k);
      for (std::size_t i = k + 1; i < n; ++i) {
        const T factor = temp(i, k) / pivot;
        if (factor == T{0}) continue;
        auto ri = temp.row(i);
        for (std::size_t j = k + 1; j < n; ++j) {
          ri[j] -= factor * rk[j];
        }
      }
    }
    return det_val;
  }

  // --- Inversão Nodal in situ (Shipley-Coleman com Pivotação Parcial) ---
  void invert_inplace() {
    if (!is_square()) {
      throw std::invalid_argument("A matriz não é quadrada");
    }

    const std::size_t n = rows_;
    if (n == 0) return;
    if (n == 1) {
      using real_t = real_type_t<T>;
      using std::abs;
      if (abs((*this)(0, 0)) < std::numeric_limits<real_t>::epsilon()) {
        throw std::runtime_error("Matriz singular: pivô nulo");
      }
      (*this)(0, 0) = T{1} / (*this)(0, 0);
      return;
    }

    using real_t = real_type_t<T>;
    using std::abs;
    const real_t eps = std::numeric_limits<real_t>::epsilon() * 100;

    std::vector<std::pair<std::size_t, std::size_t>> row_swaps;
    row_swaps.reserve(n);
    std::vector<T> col(n);

    for (std::size_t k = 0; k < n; ++k) {
      // 1. Pivotação parcial: elemento de maior magnitude na coluna k (de k a n-1)
      std::size_t pivot_row = k;
      real_t max_val = abs((*this)(k, k));
      for (std::size_t i = k + 1; i < n; ++i) {
        real_t val = abs((*this)(i, k));
        if (val > max_val) {
          max_val = val;
          pivot_row = i;
        }
      }

      if (max_val < eps) {
        throw std::runtime_error("Matriz singular ou pivô nulo na inversão");
      }

      if (pivot_row != k) {
        swap_rows(k, pivot_row);
        row_swaps.emplace_back(k, pivot_row);
      }

      // 2. Passo nodal de Shipley-Coleman
      const T pivot = (*this)(k, k);
      const T p = T{1} / pivot;

      const auto rk = row(k);
      for (std::size_t i = 0; i < n; ++i) {
        col[i] = (*this)(i, k) * p;
      }
      col[k] = T{0};

      for (std::size_t i = 0; i < n; ++i) {
        if (i == k) continue;
        const T f = col[i];
        if (f == T{0}) continue;

        auto ri = row(i);
        std::transform(ri.begin(), ri.end(), rk.begin(), ri.begin(),
                       [f](const T& y, const T& x) { return y - f * x; });
        ri[k] = f;
      }

      auto rkw = row(k);
      std::ranges::transform(rkw, rkw.begin(), [p](const T& x) { return x * p; });
      rkw[k] = -p;
    }

    // 3. Negação in situ de toda a matriz
    std::ranges::transform(data_, data_.begin(), std::negate<>{});

    // 4. Desfazer as permutações de linhas através de permutações nas colunas em ordem inversa
    for (const auto& [r1, r2] : std::ranges::reverse_view(row_swaps)) {
      swap_cols(r1, r2);
    }
  }

  [[nodiscard]] Matrix<T> inverse() const {
    Matrix<T> res = *this;
    res.invert_inplace();
    return res;
  }

  [[nodiscard]] Matrix<T> inv() const { return inverse(); }

  // --- Inversão por Adjunta / Matriz de Cofatores (Regra de Cramer) ---
  [[nodiscard]] Matrix<T> inv_adjugate() const {
    if (!is_square()) {
      throw std::invalid_argument("A matriz não é quadrada");
    }
    const std::size_t n = rows_;
    if (n == 0) return *this;

    const T det_val = det();
    using real_t = real_type_t<T>;
    using std::abs;
    if (abs(det_val) < std::numeric_limits<real_t>::epsilon() * 100) {
      throw std::runtime_error("Matriz Singular impossível inverter");
    }

    if (n == 1) {
      return Matrix<T>{{T{1} / det_val}};
    }

    Matrix<T> res(n, n);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        T cof_det = submatrix(i, j).det();
        T sign = ((i + j) % 2 == 0) ? T{1} : T{-1};
        // A matriz adjunta é a transposta da matriz dos cofatores: (adj A)_ji = C_ij
        res(j, i) = (sign * cof_det) / det_val;
      }
    }
    return res;
  }

  [[nodiscard]] Matrix<T> Inv() const { return inv_adjugate(); }

 private:
  std::size_t rows_{0};
  std::size_t cols_{0};
  std::vector<T> data_{};
};

// --- Operadores Aritméticos Livres ---

template <Scalar T>
[[nodiscard]] Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
  Matrix<T> res = lhs;
  res += rhs;
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs) {
  Matrix<T> res = lhs;
  res -= rhs;
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator-(const Matrix<T>& m) {
  Matrix<T> res = m;
  std::ranges::transform(res, res.begin(), std::negate<>{});
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B) {
  if (A.cols() != B.rows()) {
    throw std::invalid_argument("Dimensões incompatíveis para multiplicação matricial");
  }

  Matrix<T> C(A.rows(), B.cols(), T{0});
  for (std::size_t i = 0; i < A.rows(); ++i) {
    const auto rA = A.row(i);
    auto rC = C.row(i);
    for (std::size_t k = 0; k < A.cols(); ++k) {
      const T aik = rA[k];
      if (aik == T{0}) continue;
      const auto rB = B.row(k);
      for (std::size_t j = 0; j < B.cols(); ++j) {
        rC[j] += aik * rB[j];
      }
    }
  }
  return C;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator+(const Matrix<T>& m, const std::type_identity_t<T>& s) {
  Matrix<T> res = m;
  res += s;
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator+(const std::type_identity_t<T>& s, const Matrix<T>& m) {
  Matrix<T> res = m;
  res += s;
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator-(const Matrix<T>& m, const std::type_identity_t<T>& s) {
  Matrix<T> res = m;
  res -= s;
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator-(const std::type_identity_t<T>& s, const Matrix<T>& m) {
  Matrix<T> res(m.rows(), m.cols());
  for (std::size_t i = 0; i < m.rows(); ++i) {
    for (std::size_t j = 0; j < m.cols(); ++j) {
      res(i, j) = s - m(i, j);
    }
  }
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator*(const Matrix<T>& m, const std::type_identity_t<T>& s) {
  Matrix<T> res = m;
  res *= s;
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator*(const std::type_identity_t<T>& s, const Matrix<T>& m) {
  Matrix<T> res = m;
  res *= s;
  return res;
}

template <Scalar T>
[[nodiscard]] Matrix<T> operator/(const Matrix<T>& m, const std::type_identity_t<T>& s) {
  Matrix<T> res = m;
  res /= s;
  return res;
}

// --- Funções Livres Utilitárias ---

template <Scalar T>
[[nodiscard]] inline T det(const Matrix<T>& m) {
  return m.det();
}

template <Scalar T>
[[nodiscard]] inline Matrix<T> inverse(const Matrix<T>& m) {
  return m.inverse();
}

template <Scalar T>
[[nodiscard]] inline Matrix<T> transpose(const Matrix<T>& m) {
  return m.transpose();
}

template <Scalar T>
[[nodiscard]] inline Matrix<T> identity(std::size_t n) {
  return Matrix<T>::identity(n);
}

template <Scalar T>
[[nodiscard]] inline Matrix<T> eye(std::size_t n) {
  return Matrix<T>::identity(n);
}

// --- Formatação para std::ostream ---

template <Scalar T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) {
  if (m.rows() == 0 || m.cols() == 0) {
    return os << "[]\n";
  }

  std::vector<std::vector<std::string>> str_cells(m.rows(), std::vector<std::string>(m.cols()));
  std::vector<std::size_t> col_widths(m.cols(), 0);

  for (std::size_t i = 0; i < m.rows(); ++i) {
    for (std::size_t j = 0; j < m.cols(); ++j) {
      std::ostringstream ss;
      ss.copyfmt(os);
      auto val = m(i, j);
      if (std::abs(val) < 1e-14) {
        val = T{0};
      }
      ss << val;
      str_cells[i][j] = ss.str();
      if (str_cells[i][j].length() > col_widths[j]) {
        col_widths[j] = str_cells[i][j].length();
      }
    }
  }

  for (std::size_t i = 0; i < m.rows(); ++i) {
    os << "[ ";
    for (std::size_t j = 0; j < m.cols(); ++j) {
      os << std::setw(col_widths[j]) << str_cells[i][j];
      if (j + 1 < m.cols()) {
        os << "  ";
      }
    }
    os << " ]\n";
  }
  return os;
}

// Alias em português no namespace mat
template <Scalar T>
using matriz = Matrix<T>;

}  // namespace mat

// Exportação do template matriz para compatibilidade direta com código existente
template <mat::Scalar T>
using matriz = mat::Matrix<T>;

// Estrutura para armazenar o resultado da decomposição LUP (P * A = L * U)
template <mat::Scalar T>
struct LUResult {
    mat::Matrix<T> L;            // Matriz triangular inferior unitária (diag = 1)
    mat::Matrix<T> U;            // Matriz triangular superior
    mat::Matrix<T> P;            // Matriz de permutação
    std::vector<std::size_t> p;  // Vetor de permutações de linha
};

// Função para decomposição LU com pivotação parcial (PA = LU)
template <mat::Scalar T>
LUResult<T> decomposicao_lu(const mat::Matrix<T>& A) {
    if (!A.is_square()) {
        throw std::invalid_argument("A decomposição LU requer uma matriz quadrada.");
    }
    const std::size_t n = A.rows();
    mat::Matrix<T> L = mat::Matrix<T>::identity(n);
    mat::Matrix<T> U = A;
    std::vector<std::size_t> p(n);
    for (std::size_t i = 0; i < n; ++i) p[i] = i;

    using real_t = mat::real_type_t<T>;
    using std::abs;
    const real_t eps = std::numeric_limits<real_t>::epsilon() * 100;

    for (std::size_t k = 0; k < n; ++k) {
        // Encontra o pivô de maior módulo na coluna k
        std::size_t pivot_row = k;
        real_t max_val = abs(U(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            real_t val = abs(U(i, k));
            if (val > max_val) {
                max_val = val;
                pivot_row = i;
            }
        }

        if (max_val < eps) {
            throw std::runtime_error("Matriz é singular ou quase singular para decomposição LU.");
        }

        // Permuta linhas em U, no vetor de pivôs p e nas colunas já computadas de L
        if (pivot_row != k) {
            U.swap_rows(k, pivot_row);
            std::swap(p[k], p[pivot_row]);
            for (std::size_t j = 0; j < k; ++j) {
                std::swap(L(k, j), L(pivot_row, j));
            }
        }

        // Eliminação Gaussiana
        const T pivot = U(k, k);
        for (std::size_t i = k + 1; i < n; ++i) {
            T factor = U(i, k) / pivot;
            L(i, k) = factor;
            U(i, k) = T{0};
            for (std::size_t j = k + 1; j < n; ++j) {
                U(i, j) -= factor * U(k, j);
            }
        }
    }

    // Monta matriz de permutação P
    mat::Matrix<T> P(n, n, T{0});
    for (std::size_t i = 0; i < n; ++i) {
        P(i, p[i]) = T{1};
    }

    return {L, U, P, p};
}

// Inversão de matriz utilizando a decomposição LU:
// A * X = I  <=>  P * A * X = P * I  <=>  L * U * X = P
// 1) Resolve L * Y = P (substituição progressiva / forward substitution)
// 2) Resolve U * X = Y (substituição regressiva / back substitution)
template <mat::Scalar T>
mat::Matrix<T> inv_lu(const mat::Matrix<T>& A) {
    if (!A.is_square()) {
        throw std::invalid_argument("A matriz precisa ser quadrada para inversão.");
    }
    const std::size_t n = A.rows();
    auto [L, U, P, p] = decomposicao_lu(A);

    // 1. Substituição progressiva: L * Y = P
    mat::Matrix<T> Y(n, n, T{0});
    for (std::size_t col = 0; col < n; ++col) {
        for (std::size_t i = 0; i < n; ++i) {
            T sum = P(i, col);
            for (std::size_t j = 0; j < i; ++j) {
                sum -= L(i, j) * Y(j, col);
            }
            Y(i, col) = sum; // L(i, i) == 1
        }
    }

    // 2. Substituição regressiva: U * X = Y
    mat::Matrix<T> X(n, n, T{0});
    for (std::size_t col = 0; col < n; ++col) {
        for (std::size_t step = 0; step < n; ++step) {
            std::size_t i = n - 1 - step;
            T sum = Y(i, col);
            for (std::size_t j = i + 1; j < n; ++j) {
                sum -= U(i, j) * X(j, col);
            }
            X(i, col) = sum / U(i, i);
        }
    }

    // 3. Refinamento Iterativo de Newton-Schultz: eleva a precisão ao limite da máquina
    mat::Matrix<T> I = mat::Matrix<T>::identity(n);
    mat::Matrix<T> E = I - A * X;
    X += X * E;

    return X.chop(1e-14);
}


