# Solution of Linear Equations (Iterative Methods)

In numerical analysis, solving systems of linear equations—typically written as

\[
A\mathbf{x} = \mathbf{b}
\]

is a fundamental problem. While **direct methods** such as *Cramer’s Rule* or *Gaussian Elimination* aim to compute an exact solution in a finite number of steps, they become computationally expensive and memory-intensive for large-scale systems.

This limitation motivates the use of **iterative methods**, which start from an initial guess and progressively refine the solution until the error becomes acceptably small.

---

## 1. Iterative Methods: The Art of Successive Refinement

### Why “Iterative”?

The term *iterative* comes from the Latin *iterare*, meaning *to repeat*. Instead of solving the system in a single pass, iterative algorithms generate a sequence of approximations:

\[
\{\mathbf{x}^{(0)}, \mathbf{x}^{(1)}, \mathbf{x}^{(2)}, \ldots\}
\]

Each new approximation is computed using information from the previous one, gradually reducing the error.

> **Intuition:**  
> Iterative methods are similar to tuning a guitar string—pluck, listen, adjust slightly, and repeat until the desired pitch is reached.

Iterative methods are particularly effective for:
- Large and sparse systems
- Problems where approximate solutions are sufficient
- Parallel and memory-efficient computations

---

## (i) Jacobi Iterative Method

### Theory: Simultaneous Displacement

The **Jacobi method** is one of the simplest iterative techniques. Each equation of the system is solved for its corresponding variable \(x_i\), assuming all other variables are known from the **previous iteration**.

The defining characteristic of the Jacobi method is that **all updates are computed simultaneously** using only old values:

\[
x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{\substack{j=1 \\ j \ne i}}^{n} a_{ij} x_j^{(k)} \right)
\]

Because no updated value is reused within the same iteration, the Jacobi method is easy to implement and highly suitable for **parallel computation**.

---

### Convergence Requirement

For guaranteed convergence, the coefficient matrix \(A\) should be **diagonally dominant**:

\[
|a_{ii}| > \sum_{j \ne i} |a_{ij}| \quad \text{for all } i
\]

---

### Algorithm

1. **Rearrangement**  
   Rewrite each equation so that \(x_i\) appears on the left-hand side.

2. **Initial Guess**  
   Choose an initial vector \(\mathbf{x}^{(0)}\) (commonly all zeros).

3. **Iteration**  
   Compute each \(x_i^{(k+1)}\) using only values from iteration \(k\).

4. **Stopping Criterion**  
   Stop when:
   \[
   \|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\| < \varepsilon
   \]
   where \(\varepsilon\) is a prescribed tolerance.

---

### Pseudocode

```plaintext
Input: Matrix A, Vector b, tolerance e, max_iterations N

Initialize:
    x_old = [0, 0, ..., 0]
    x_new = [0, 0, ..., 0]

For k from 1 to N:
    For i from 1 to rows(A):
        sum = 0
        For j from 1 to columns(A):
            if i ≠ j:
                sum = sum + A[i][j] * x_old[j]

        x_new[i] = (b[i] - sum) / A[i][i]

    If distance(x_new, x_old) < e:
        Print "Converged"
        Break

    x_old = x_new   // prepare for next iteration
