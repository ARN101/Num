## 2. Solution of Linear Equations: Iterative Methods

### Overview
Iterative methods are algorithms for solving systems of linear equations ($Ax = b$) by generating a sequence of approximate solutions that converge to the exact solution. Unlike direct methods (e.g., Gaussian Elimination), which aim to solve the system in a finite number of steps, iterative methods start with an initial guess and successively improve it. They are particularly efficient for **large, sparse matrices** where direct methods would be computationally expensive or memory-intensive.

Two fundamental iterative techniques are the **Jacobi Method** and the **Gauss-Seidel Method**. Both require the matrix $A$ to be **diagonally dominant** to guarantee convergence.

### A. Jacobi Iterative Method

The Jacobi method is a method of **simultaneous displacement**. To find the $k+1$ approximation for the variable $x_i$, the method uses **only** the values from the previous iteration $k$.

**Formula:**


$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1, j \ne i}^{n} a_{ij} x_j^{(k)} \right), \quad i = 1, \dots, n$$


#### Algorithm
1.  **Input:** Matrix $A$, Vector $b$, Initial guess $x^{(0)}$, Tolerance $\epsilon$, Max Iterations $N$.
2.  **Initialize:** Create a temporary array `x_new` to store updated values.
3.  **Loop ($k = 0$ to $N$):**
    * For each variable $i$ from 1 to $n$:
        * Calculate sum $\sigma = \sum_{j \ne i} a_{ij} x_j^{(k)}$ using values from the **current** $x$ array.
        * Compute $x\_new[i] = (b[i] - \sigma) / a[i][i]$.
    * **Calculate Error:** $Error = ||x\_new - x||$.
    * **Update:** Set $x = x\_new$.
    * **Stop Condition:** If $Error < \epsilon$, break loop.
4.  **Output:** Solution vector $x$.

#### Pseudocode
```text
FUNCTION Jacobi(A, b, x0, tol, max_iter)
    n = length(b)
    x = x0
    
    FOR k = 1 TO max_iter DO
        x_new = copy(x) // Important: Store updates in a separate array
        
        FOR i = 0 TO n-1 DO
            sum = 0
            FOR j = 0 TO n-1 DO
                IF j != i THEN
                    sum = sum + A[i][j] * x[j] // Uses OLD x values
                END IF
            END FOR
            
            x_new[i] = (b[i] - sum) / A[i][i]
        END FOR
        
        // Check for convergence
        max_diff = 0
        FOR i = 0 TO n-1 DO
            max_diff = MAX(max_diff, ABS(x_new[i] - x[i]))
        END FOR
        
        x = x_new // Update x for the next iteration
        
        IF (max_diff < tol) THEN
            RETURN x
        END IF
    END FOR
    
    PRINT "Max iterations reached"
    RETURN x
END FUNCTION
```

### B. Gauss-Seidel Iterative Method

#### Mathematical Theory
The Gauss-Seidel method is a method of **successive displacement**. It updates the solution vector in place. When calculating $x_i^{(k+1)}$, it uses the **newly computed** values for $x_1, \dots, x_{i-1}$ from the *current* iteration ($k+1$), and old values for the rest ($x_{i+1}, \dots, x_n$). This typically results in faster convergence than Jacobi.

**Formula:**
$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k)} \right), \quad i = 1, \dots, n$$

#### Algorithm
1.  **Input:** Matrix $A$, Vector $b$, Initial guess $x$, Tolerance $\epsilon$.
2.  **Loop:**
    * For each variable $i$ from 1 to $n$:
        * Calculate sum $\sigma$ using the **most recent** values of $x$ available in memory.
        * Compute new value: $val = (b[i] - \sigma) / a[i][i]$.
        * **Update Immediately:** $x[i] = val$ (overwrite the old value).
    * **Check Convergence:** If the maximum change in any $x[i]$ is less than $\epsilon$, stop.
3.  **Output:** Solution vector $x$.

### Real-World Applications
1.  **Power System Load Flow:** In electrical engineering, the Gauss-Seidel method is a classic approach for solving the non-linear power flow equations to determine voltage magnitudes and angles at different buses in a power grid.
2.  **Computer Graphics (Fluid Simulation):** Jacobi iterations are often used in solving the pressure projection step in the Navier-Stokes equations for real-time fluid dynamics (e.g., in video games) because the independence of updates makes it highly parallelizable on GPUs.
3.  **Thermal Analysis:** Solving the steady-state heat distribution in a solid object involves discretizing Laplace’s equation ($\nabla^2 T = 0$) into a large system of linear equations, which are efficiently solved using these iterative solvers.

### Further Reading
* **Comparison of Methods:** [Wolfram MathWorld - Jacobi Method](https://mathworld.wolfram.com/JacobiMethod.html) | [Gauss-Seidel Method](https://mathworld.wolfram.com/Gauss-SeidelMethod.html)
* **Deep Dive into Convergence:** [MIT OpenCourseWare - Matrix Methods in Data Analysis](https://ocw.mit.edu/courses/18-065-matrix-methods-in-data-analysis-signal-processing-and-machine-learning-spring-2018/)
* **Numerical Recipes:** [Iterative Solution of Linear Algebraic Equations](http://numerical.recipes/book/book.html)
