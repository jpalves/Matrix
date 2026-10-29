# Guia Teórico e Prático: Decomposição LU, Inversão e SVD

Este documento explica em detalhes os fundamentos matemáticos, algoritmos numéricos e decisões de implementação dos módulos de **Decomposição LU**, **Inversão Matricial** e **Decomposição em Valores Singulares (SVD)** implementados no projeto (`mat.hpp`, `svd.hpp`, `complexo.hpp`, `inv.cpp`).

---

## Sumário
1. [Visão Geral da Arquitetura](#1-visão-geral-da-arquitetura)
2. [Decomposição LU com Pivotação Parcial ($PA = LU$)](#2-decomposição-lu-com-pivotação-parcial-pa--lu)
   - [O Problema dos Pivôs Nulos](#o-problema-dos-pivôs-nulos)
   - [Algoritmo de Doolittle com Pivotação Parcial](#algoritmo-de-doolittle-com-pivotação-parcial)
   - [Convenção da Matriz de Permutação ($PA = LU$ vs $A = PLU$)](#convenção-da-matriz-de-permutação)
   - [Inversão de Matrizes via Decomposição LU](#inversão-de-matrizes-via-decomposição-lu)
   - [Pivotação de Linhas vs Pivotação de Colunas](#pivotação-de-linhas-vs-pivotação-de-colunas)
3. [Decomposição em Valores Singulares (SVD)](#3-decomposição-em-valores-singulares-svd)
   - [Teorema Fundamental](#teorema-fundamental)
   - [O Método One-Sided Jacobi (Hestenes)](#o-método-one-sided-jacobi-hestenes)
   - [Formulação para Matrizes Reais](#formulação-para-matrizes-reais)
   - [Formulação para Matrizes Complexas (Rotações Unitárias)](#formulação-para-matrizes-complexas)
   - [Tratamento de Matrizes Retangulares ($m < n$ e $m > n$)](#matrizes-retangulares)
4. [Inversão e Pseudo-Inversa de Moore-Penrose via SVD](#4-inversão-e-pseudo-inversa-de-moore-penrose-via-svd)
5. [Controle de Estabilidade e Redução de Erro Numérico](#5-controle-de-estabilidade-e-redução-de-erro-numérico)
   - [A Origem do Ruído Numérico ($\approx 10^{-16}$)](#a-origem-do-ruído-numérico)
   - [Refinamento Iterativo de Newton-Schultz](#refinamento-iterativo-de-newton-schultz)
   - [Cálculo em Precisão Estendida (`long double`)](#cálculo-em-precisão-estendida)
   - [Filtragem de Ruído (`chop` / `limpa`)](#filtragem-de-ruído-chop)
6. [Exemplos de Código (C++ e Python)](#6-exemplos-de-código)

---

## 1. Visão Geral da Arquitetura

O projeto organiza a computação matricial em três cabeçalhos modulares em C++20:
- **`mat.hpp`**: Classe de matriz densa `mat::Matrix<T>` (ou `matriz<T>`) com suporte genérico a escalares reais (`float`, `double`) e complexos (`std::complex<T>`, `complexo<T>`). Inclui transposição, transposição conjugada (Hermitiana), determinante, e inversão por LU.
- **`complexo.hpp`**: Classe que estende `std::complex<T>` com sobrecargas aritméticas mistas (`double * complexo`, `complexo / double`, etc.) e formatação.
- **`svd.hpp`**: Implementação autônoma do algoritmo One-Sided Jacobi SVD para matrizes reais e complexas, pseudo-inversa `pinv` e utilitários como `diag`.

---

## 2. Decomposição LU com Pivotação Parcial ($PA = LU$)

### O Problema dos Pivôs Nulos
A decomposição direta $A = L U$ consiste em fatorar uma matriz quadrada $A$ em:
- $L$: triangular inferior com $1$s na diagonal principal ($L_{i,i} = 1$).
- $U$: triangular superior.

Porém, na eliminação Gaussiana clássica, os multiplicadores são calculados por:
$$m_{i, k} = \frac{U_{i, k}}{U_{k, k}}$$
Se em qualquer etapa o pivô $U_{k, k}$ for nulo ou próximo de zero (como na matriz do exemplo, cuja linha 1 é `[0, 0, 1]`), ocorre **divisão por zero** ou instabilidade numérica catastrófica.

### Algoritmo de Doolittle com Pivotação Parcial
Para resolver esse problema, introduz-se a matriz de permutação $P$, de forma que:
$$P \cdot A = L \cdot U$$

Em cada coluna $k \in [0, n-1]$:
1. **Busca do Pivô**: Seleciona-se a linha $p \ge k$ com maior valor absoluto na coluna $k$:
   $$p = \arg\max_{i \ge k} |U_{i, k}|$$
2. **Troca de Linhas**:
   - Troca as linhas $k$ e $p$ da matriz $U$.
   - Troca as entradas correspondentes no vetor de permutações $p$.
   - Troca os multiplicadores já computados em $L$ nas colunas anteriores ($j < k$).
3. **Eliminação**:
   Para cada linha $i > k$:
   $$\text{fator} = \frac{U_{i, k}}{U_{k, k}}, \quad L_{i, k} = \text{fator}, \quad U_{i, k} = 0$$
   $$U_{i, j} \leftarrow U_{i, j} - \text{fator} \cdot U_{k, j} \quad (\forall j > k)$$

### Convenção da Matriz de Permutação
Existem duas convenções principais na literatura:
1. **Convenção $P A = L U$ (MATLAB, nosso código C++):**
   $P$ multiplica $A$ à esquerda, permutando as **linhas** de $A$ para que ela seja fatorável em $L U$.
2. **Convenção $A = P_{\text{scipy}} L U$ (SciPy / LAPACK):**
   Aqui $P_{\text{scipy}} = P^{-1} = P^T$.
   Como a transposta de uma matriz de permutação permuta linhas com colunas, **as colunas da matriz $P$ do SciPy aparecem trocadas em relação ao $P$ de $P A = L U$**.

Ambas as formas são matematicamente equivalentes através de $P_{\text{scipy}} = P^T$.

### Inversão de Matrizes via Decomposição LU
Para inverter uma matriz $A$ de ordem $n$, queremos resolver o sistema matricial para a matriz de incógnitas $X$:
$$A \cdot X = I$$

Multiplicando à esquerda por $P$:
$$P A X = P I \implies L \cdot U \cdot X = P$$

O problema decompõe-se em duas etapas triangulares muito rápidas ($\mathcal{O}(n^2)$ por coluna):
1. **Substituição Progressiva (Forward Substitution):**
   Resolve $L \cdot Y = P$ coluna por coluna:
   $$Y_{i, \text{col}} = P_{i, \text{col}} - \sum_{j=0}^{i-1} L_{i, j} \cdot Y_{j, \text{col}} \quad (\text{pois } L_{i, i} = 1)$$
2. **Substituição Regressiva (Back Substitution):**
   Resolve $U \cdot X = Y$ coluna por coluna, de baixo para cima:
   $$X_{i, \text{col}} = \frac{1}{U_{i, i}} \left( Y_{i, \text{col}} - \sum_{j=i+1}^{n-1} U_{i, j} \cdot X_{j, \text{col}} \right)$$

Ao final, a matriz $X$ é exatamente $A^{-1}$.

### Pivotação de Linhas vs Pivotação de Colunas
Considere a matriz do nosso teste:
$$A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 7 & 1 & 2 \end{pmatrix}$$

- **Pivotação de Linhas ($P A = L U$):** A eliminação busca pivôs variando as linhas na coluna atual. Ela troca a linha 0 com a linha 2, e depois a linha 1 com a 2.
- **Pivotação de Colunas ($A Q = L U$):** Se trocarmos a coluna 1 com a coluna 2 de $A$ via matriz de permutação $Q$:
  $$A \cdot Q = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 7 & 2 & 1 \end{pmatrix}$$
  Esta matriz **já é uma matriz triangular inferior perfeita ($L$)**! Portanto:
  $$A Q = L \cdot I \implies A = L \cdot Q^T$$
  Ou seja, a matriz $A$ é intrinsecamente uma triangular inferior com as colunas 1 e 2 trocadas.

---

## 3. Decomposição em Valores Singulares (SVD)

### Teorema Fundamental
Qualquer matriz $A \in \mathbb{C}^{m \times n}$ pode ser decomposta como:
$$A = U \cdot \Sigma \cdot V^H$$
- $U \in \mathbb{C}^{m \times m}$: matriz unitária ($U^H U = I$). Suas colunas são os vetores singulares à esquerda.
- $\Sigma \in \mathbb{R}^{m \times n}$: matriz diagonal retangular contendo os valores singulares $\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_r \ge 0$.
- $V \in \mathbb{C}^{n \times n}$: matriz unitária ($V^H V = I$), onde $V^H$ é a transposta conjugada (Hermitiana). Suas colunas são os vetores singulares à direita.

### O Método One-Sided Jacobi (Hestenes)
Diferente de métodos baseados em bidiagonalização QR (que exigem passos complexos e bibliotecas pesadas), o **Método de Hestenes (One-Sided Jacobi)** é auto-contido e possui a maior precisão relativa comprovada para valores singulares de diferentes ordens de grandeza (Teorema de Demmel-Kahan).

**Ideia central:** Em vez de diagonalizar $A$ pelos dois lados simultaneamente, aplicamos rotações unitárias $J$ à direita de $A$:
$$B = A \cdot V$$
Iteramos até que todas as colunas de $B$ sejam **mutuamente ortogonais**:
$$\langle B_{*, j}, B_{*, k} \rangle = 0 \quad (\forall j \neq k)$$

Quando a convergência é atingida:
- O comprimento da coluna $j$ de $B$ é o valor singular: $\sigma_j = \|B_{*, j}\|_2$.
- Normalizando cada coluna não-nula, obtemos $U$: $U_{*, j} = \frac{B_{*, j}}{\sigma_j}$.
- A matriz acumulada de rotações é $V$, e $A = U \Sigma V^H$.

---

### Formulação para Matrizes Reais
Para cada par de colunas $(j, k)$ com $j < k$:
1. Computamos os produtos internos:
   $$\alpha = \|col_j\|^2, \quad \beta = \|col_k\|^2, \quad \gamma = col_j \cdot col_k$$
2. Queremos uma rotação plana no plano $(j, k)$:
   $$\begin{pmatrix} col_j' & col_k' \end{pmatrix} = \begin{pmatrix} col_j & col_k \end{pmatrix} \begin{pmatrix} c & -s \\ s & c \end{pmatrix}$$
   tal que $col_j' \cdot col_k' = 0$.
3. Desenvolvendo a equação:
   $$(c^2 - s^2) \gamma + c s (\beta - \alpha) = 0 \implies \tan(2\theta) = \frac{2\gamma}{\alpha - \beta}$$
4. Definindo $\zeta = \frac{\alpha - \beta}{2\gamma}$, a raiz estável para $t = \tan\theta$ é:
   $$t = \frac{\text{sinal}(\zeta)}{|\zeta| + \sqrt{1 + \zeta^2}}$$
   $$c = \frac{1}{\sqrt{1 + t^2}}, \quad s = t \cdot c$$
5. Atualizam-se as colunas $j$ e $k$ de $B$ e da matriz $V$.

---

### Formulação para Matrizes Complexas

Para matrizes com elementos complexos, três diferenças fundamentais devem ser respeitadas:

#### 1. Produto Interno Hermitiano
O produto interno entre as colunas é complexo:
$$\gamma = col_j^H col_k = \sum_{i=0}^{m-1} \overline{B_{i, j}} \cdot B_{i, k} \in \mathbb{C}$$
enquanto as normas $\alpha = \sum |B_{i, j}|^2$ e $\beta = \sum |B_{i, k}|^2$ são reais.

#### 2. Rotação Unitária Complexa ($2 \times 2$)
A matriz de transformação deve ser unitária no plano $(j, k)$:
$$J = \begin{pmatrix} c & s \\ -\overline{s} & c \end{pmatrix}$$
com $c \in \mathbb{R}$ e $s \in \mathbb{C}$, garantindo $c^2 + |s|^2 = 1$ e $J^H J = I$.

As colunas são atualizadas por:
$$col_j' = c \cdot col_j + s \cdot col_k$$
$$col_k' = -\overline{s} \cdot col_j + c \cdot col_k$$

#### 3. Cancelamento da Fase Complexa
Impondo a ortogonalidade hermitiana $(col_j')^H col_k' = 0$:
$$(c col_j^H + \overline{s} col_k^H)(-\overline{s} col_j + c col_k) = 0$$
$$-c \overline{s} \alpha + c^2 \gamma - (\overline{s})^2 \overline{\gamma} + c \overline{s} \beta = 0$$

Ao escolhermos o parâmetro complexo com a fase oposta à de $\gamma$:
$$s = s_0 \cdot \frac{\overline{\gamma}}{|\gamma|}, \quad \text{onde } s_0 \in \mathbb{R}$$
temos $\overline{s} = s_0 \frac{\gamma}{|\gamma|}$. Substituindo na equação, a fase $\frac{\gamma}{|\gamma|}$ fatora completamente:
$$\frac{\gamma}{|\gamma|} \left[ (c^2 - s_0^2) |\gamma| - c s_0 (\alpha - \beta) \right] = 0$$

A equação restante para $c$ e $s_0$ torna-se idêntica à do caso real, substituindo $\gamma$ por seu módulo $|\gamma|$:
$$\zeta = \frac{\alpha - \beta}{2 |\gamma|}, \quad t = \frac{\text{sinal}(\zeta)}{|\zeta| + \sqrt{1 + \zeta^2}}$$
$$c = \frac{1}{\sqrt{1 + t^2}}, \quad s_0 = t \cdot c$$
$$s = s_0 \cdot \frac{\overline{\gamma}}{|\gamma|}$$

Isso garante convergência para qualquer matriz complexa mantendo $U$ e $V$ estritamente unitárias.

---

### Matrizes Retangulares
- **Se $m \ge n$:** O algoritmo produz diretamente $U_{m \times n}$, $\Sigma_{n \times n}$ e $V^H_{n \times n}$. Se necessária a base completa $U_{m \times m}$, complementam-se as $m - n$ colunas restantes por ortogonalização de Gram-Schmidt.
- **Se $m < n$:** Aplica-se o SVD à matriz adjunta $A^H$:
  $$A^H = U_t \Sigma_t V_t^H \implies A = (U_t \Sigma_t V_t^H)^H = V_t \Sigma_t^T U_t^H$$
  Logo, $U = V_t$ e $V^H = U_t^H$.

---

## 4. Inversão e Pseudo-Inversa de Moore-Penrose via SVD

A partir de $A = U \Sigma V^H$, a pseudo-inversa $A^+$ é obtida invertendo os valores singulares não-nulos:
$$A^+ = V \cdot \Sigma^+ \cdot U^H$$
onde:
$$\Sigma^+_{i, i} = \begin{cases} \frac{1}{\sigma_i}, & \text{se } \sigma_i > \text{tol} \cdot \sigma_1 \\ 0, & \text{caso contrário} \end{cases}$$

Se $A$ for quadrada e inversível, $A^+ = A^{-1}$. Para matrizes singulares ou retangulares, $A^+$ fornece a solução de mínimos quadrados com norma mínima.

---

## 5. Controle de Estabilidade e Redução de Erro Numérico

### A Origem do Ruído Numérico
Em aritmética de ponto flutuante IEEE 754 de precisão dupla (`double`), o menor incremento relativo é:
$$\varepsilon_{\text{máquina}} \approx 2.22 \times 10^{-16}$$

Quando produtos internos ou subtrações teóricas deveriam resultar em zero, pequenas diferenças de arredondamento nos bits menos significativos produzem resíduos da ordem de $\pm 5.55 \times 10^{-17}$ ou $\pm 6.93 \times 10^{-18}$.

### Refinamento Iterativo de Newton-Schultz
Para eliminar esse resíduo e levar a inversa à máxima exatidão teórica da máquina, implementamos o **método de Newton-Schultz**:
Dada uma aproximação $X_0 \approx A^{-1}$, define-se a matriz de resíduo:
$$E = I - A \cdot X_0$$
A atualização iterativa é:
$$X_1 = X_0 + X_0 \cdot E = X_0 (2I - A X_0)$$

O novo resíduo é:
$$I - A X_1 = I - A X_0 (I + E) = I - (I - E)(I + E) = E^2$$
Como $\|E\| \approx 10^{-16}$, o novo resíduo passa a ser:
$$\|E^2\| \approx (10^{-16})^2 = 10^{-32}$$
Isso reduz o resíduo abaixo da resolução do `double`, alcançando exatidão plena.

### Cálculo em Precisão Estendida
No cabeçalho `svd.hpp`, as acumulações de $\alpha, \beta, \gamma$ e os cálculos trigonométricos são processados em `long double` (80 bits no x86_64, com $\varepsilon \approx 1.08 \times 10^{-19}$), prevenindo o cancelamento numérico catastrófico.

### Filtragem de Ruído (`chop` / `limpa`)
A classe `mat::Matrix` disponibiliza o método:
```cpp
M.chop(1e-14); // zera qualquer componente real ou imag com módulo < 1e-14
```
Além disso, o operador de visualização `operator<<` formata automaticamente resíduos com $|x| < 10^{-14}$ como `0`.

---

## 6. Exemplos de Código

### Exemplo em C++ (`inv.cpp`)
```cpp
#include <iostream>
#include "mat.hpp"
#include "complexo.hpp"
#include "svd.hpp"

using namespace std;

int main() {
    complexo<double> j = sqrt(-1);
    matriz<complexo<double>> A({{1, 0, 0},
                                {0, 0, 1},
                                {7, 1, 2}});

    // 1. Decomposição LU e Inversão
    auto lu = decomposicao_lu(A);
    matriz<complexo<double>> A_inv = inv_lu(A);
    cout << "Inversa por LU:\n" << A_inv << endl;

    // 2. SVD de Matriz Complexa (Aj = A * j)
    matriz<complexo<double>> Aj = A * j;
    auto [U, S, Vh, Vt, Sigma] = svd(Aj);

    cout << "Valores singulares de Aj:\n";
    for (auto s : S) cout << s << " ";
    cout << "\n\nReconstrução (U * Sigma * V^H):\n";
    cout << U * Sigma * Vh << endl;

    // 3. Inversa por SVD (pseudo-inversa com Newton-Schultz)
    cout << "Inversa de Aj por SVD:\n" << pinv(Aj) << endl;

    return 0;
}
```

### Exemplo Equivalente em Python (NumPy / SciPy)
```python
import numpy as np
import scipy.linalg

j = 1j
A = np.array([[1, 0, 0],
              [0, 0, 1],
              [7, 1, 2]], dtype=complex) * j

# LU com SciPy (A = P @ L @ U)
P, L, U = scipy.linalg.lu(A)
A_inv_lu = np.linalg.inv(U) @ np.linalg.inv(L) @ P.T

# SVD com NumPy (A = U @ diag(S) @ Vh)
u, s, vh = np.linalg.svd(A)
A_rec = u @ np.diag(s) @ vh
A_pinv = vh.conj().T @ np.diag(1.0 / s) @ u.conj().T

print("Valores singulares:", s)
print("Reconstrução SVD:\n", np.round(A_rec, 4))
print("Inversa SVD:\n", np.round(A_pinv, 4))
```
