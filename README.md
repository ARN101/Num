## A. Solution of Linear Equations

In numerical analysis, solving linear systems (typically written as $Ax = b$) is fundamental. While you might be familiar with "Direct Methods" like Cramer's Rule or Gaussian Elimination that attempt to find the exact solution in a finite number of steps, they can become computationally expensive for very large systems.

This is where **Iterative Methods** shine. Instead of trying to solve the problem in one go, these methods start with a guess and refine it over and over again until the error is negligible.

### 1. Iterative Methods: The Art of Successive Refinement

**Why "Iterative"?**
The term comes from the Latin *iterare* (to repeat). Unlike direct elimination, these algorithms generate a sequence of approximate solutions $\{x^{(0)}, x^{(1)}, x^{(2)}, ...\}$. Each step "iterates" on the previous one to reduce the error. Think of it like tuning a guitar: you pluck the string, check the pitch, adjust the peg slightly, and repeat until it sounds perfect.

#### (i) Jacobi Iterative Method

**Theory: Simultaneous Displacement**
The Jacobi method is the simplest iterative technique. It works by isolating the variable $x_i$ in the $i$-th equation. The unique characteristic of Jacobi is that it uses values from the **previous** iteration to calculate **all** new values. No new information is used until the next full cycle.

Because each update is independent of the others within the same step, this method is famously easy to parallelize on modern multi-core processors.

**Algorithm**
1.  **Arrangement:** Rewrite the system so that $x_1$ is on the left of equation 1, $x_2$ on the left of equation 2, etc.
    * *Note: The system must be diagonally dominant ($|a_{ii}| > \sum |a_{ij}|$) for guaranteed convergence.*
2.  **Guess:** Start with an initial guess $x^{(0)}$ (often all zeros).
3.  **Iterate:** For each variable $x_i$, compute the new value using only the old values from step $k$.
4.  **Stop:** Repeat until the difference between the new and old values is less than your allowed error tolerance ($\epsilon$).

**Pseudocode**
```text
Input: Matrix A, Vector b, tolerance e, max_iterations N
Initialize: x_old = [0, 0, ... 0]
            x_new = [0, 0, ... 0]

For k from 1 to N:
    For i from 1 to rows(A):
        sum = 0
        For j from 1 to columns(A):
            if i != j:
                sum = sum + A[i][j] * x_old[j]
        
        x_new[i] = (b[i] - sum) / A[i][i]
    
    If distance(x_new, x_old) < e:
        Print "Converged"
        Break
    
    x_old = x_new // Update for next cycle
