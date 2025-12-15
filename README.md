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
```
**Further Study**
* [Jacobi Method - Wikipedia (Comprehensive Theory)](https://en.wikipedia.org/wiki/Jacobi_method)
* [Jacobi Method Explained - GeeksforGeeks (Examples & Code)](https://www.geeksforgeeks.org/engineering-mathematics/jacobian-method/)

#### (ii) Gauss-Seidel Iterative Method

**Theory: Successive Displacement**
The Gauss-Seidel method is an optimization of the Jacobi technique. In the Jacobi method, we hold all updates in a buffer until the end of the iteration. In Gauss-Seidel, we update the variables **immediately**.

As soon as a new value for $x_1$ is calculated, it replaces the old $x_1$ in memory and is used immediately to calculate $x_2$. Because we are using the most recent information available at every step, Gauss-Seidel typically converges much faster than Jacobi (often requiring half as many iterations). However, because $x_2$ depends on the *new* $x_1$, the steps must be done in order, making parallelization difficult.

**Algorithm**
1.  **Rearrange:** Same setup as Jacobi (requires diagonal dominance).
2.  **Initialize:** Start with an initial guess vector $x$.
3.  **Iterate:** Calculate $x_i$ and immediately overwrite the value in the solution vector. Use this new value for calculations of $x_{i+1}$, $x_{i+2}$, etc.
4.  **Check Convergence:** Stop when the error falls below the tolerance level.

**Pseudocode**
```text
Input: Matrix A, Vector b, tolerance e, max_iterations N
Output: Solution vector x

Initialize: 
    x = [0, 0, ... 0] // Only one array is needed

For k from 1 to N:
    converged = true
    
    For i from 1 to rows(A):
        old_val = x[i]
        sum = 0
        For j from 1 to columns(A):
            if i != j:
                // Uses the most recent x[j] available in memory
                sum = sum + A[i][j] * x[j] 
        
        x[i] = (b[i] - sum) / A[i][i]
        
        // Check error for this specific variable
        If abs(x[i] - old_val) > e:
            converged = false
            
    If converged is true:
        Print "Converged"
        Return x
```
**Further Study**
* [Gauss-Seidel Method - GeeksforGeeks](https://www.geeksforgeeks.org/gauss-seidel-method/)
* [Iterative Methods: Jacobi vs Gauss-Seidel - Math LibreTexts](https://math.libretexts.org/Bookshelves/Linear_Algebra/Introduction_to_Matrix_Algebra_(Kaw)/01%3A_Chapters/1.08%3A_Gauss-Seidel_Method)

## D. Numerical Differentiation

Numerical differentiation is the process of calculating the derivative (rate of change) of a function using a set of discrete data points $(x_i, y_i)$ rather than an analytical formula. This is essential in real-world engineering where the underlying function is often unknown (e.g., sensor data) or too complex to differentiate manually.

Instead of taking the limit as $h \to 0$ analytically, we approximate the slope using finite steps $h$.

### 1. Equal-Interval Interpolation Method

**Theory: Differentiating the Polynomial**
When data points are spaced equally (with a constant step size $h$), we can approximate the function $f(x)$ using Newton's Interpolation formulas (Forward or Backward) and then differentiate that polynomial.

For a point $x$ near the beginning of the data set, we use **Newton's Forward Difference formula**. The first derivative is approximated as:

$$\frac{dy}{dx} = \frac{1}{h} \left[ \Delta y_0 + \frac{2p-1}{2} \Delta^2 y_0 + \frac{3p^2-6p+2}{6} \Delta^3 y_0 + \dots \right]$$

Where $p = \frac{x - x_0}{h}$.
If we are calculating the derivative exactly at a tabulated point ($x = x_0$, so $p=0$), the formula simplifies significantly to:
$$f'(x_0) \approx \frac{1}{h} \left( \Delta y_0 - \frac{1}{2}\Delta^2 y_0 + \frac{1}{3}\Delta^3 y_0 - \dots \right)$$

**Algorithm**
1.  **Check Interval:** Verify that all $x$ values have a constant difference $h$.
2.  **Difference Table:** Construct the Forward Difference Table ($\Delta$) if $x$ is near the start, or Backward Difference Table ($\nabla$) if $x$ is near the end.
3.  **Calculate $p$:** Determine the position factor $p = (x_{target} - x_0) / h$.
4.  **Apply Series:** Substitute the difference values ($\Delta y_0, \Delta^2 y_0, \dots$) and $p$ into the differentiation formula.
5.  **Scale:** Divide the result by step size $h$ to get the final derivative.

**Pseudocode**
```text
Input: Arrays x[] and y[], value xp (point to differentiate at)
Output: Derivative value dy/dx

1. Calculate step size h = x[1] - x[0]
2. Generate Forward Difference Table diff[n][n]
   For j = 1 to n-1:
       For i = 0 to n-j-1:
           diff[i][j] = diff[i+1][j-1] - diff[i][j-1]

3. Find index i such that x[i] is closest to xp
4. Calculate p = (xp - x[i]) / h
5. Initialize sum = 0

6. Apply Formula (Example for 1st derivative):
   term1 = diff[i][0]
   term2 = (2*p - 1) * diff[i][1] / 2
   term3 = (3*p*p - 6*p + 2) * diff[i][2] / 6
   ... (continue for higher orders)
   
   sum = term1 + term2 + term3 + ...

7. Result = sum / h
8. Return Result
```
**Further Study**
* [Newton’s Forward Difference Formula for Differentiation - GeeksforGeeks](https://www.geeksforgeeks.org/newtons-forward-difference-formula-for-differentiation/)
* [Numerical Differentiation - Math.OHIO.edu](https://web.math.ohio.edu.cn/~courses/math3600/Lecture13.pdf)

## F. Numerical Integration

Numerical integration, often called "numerical quadrature," is the process of calculating the approximate value of a definite integral $\int_{a}^{b} f(x) dx$. While analytical calculus finds the exact area under a curve using antiderivatives, numerical methods sum the areas of geometric shapes (like trapezoids or parabolas) fitted under the curve.

This is critical in simulations where $f(x)$ is not a simple formula but a stream of data points (e.g., calculating distance from a velocity-time log).

### 1. Simpson’s 1/3 Rule

**Theory: Parabolic Approximation**
Simpson’s 1/3 Rule improves upon the Trapezoidal Rule by approximating the function $f(x)$ not as a straight line, but as a **second-order polynomial (parabola)** connecting every three points.



Because it fits parabolas, it requires an **even number of segments** (intervals) $n$, which means you need an odd number of data points. The formula weights the boundary points and internal points differently to achieve higher accuracy.

$$I \approx \frac{h}{3} \left[ (y_0 + y_n) + 4(y_1 + y_3 + \dots + y_{n-1}) + 2(y_2 + y_4 + \dots + y_{n-2}) \right]$$

**Algorithm**
1.  **Verify Intervals:** Ensure the number of intervals $n$ is even. If $n$ is odd, this method cannot be applied directly over the whole range.
2.  **Calculate Step Size:** $h = (b - a) / n$.
3.  **Sum Extremes:** Add the first ($y_0$) and last ($y_n$) values.
4.  **Sum Odds:** Multiply the sum of ordinates at odd positions ($y_1, y_3...$) by 4.
5.  **Sum Evens:** Multiply the sum of ordinates at even positions ($y_2, y_4...$) by 2.
6.  **Calculate Total:** Sum all components and multiply by $h/3$.

**Pseudocode**
```text
Input: Function f(x), lower_limit a, upper_limit b, intervals n
Output: Integral value

If n % 2 != 0:
    Print "Error: n must be even"
    Return

h = (b - a) / n
sum = f(a) + f(b) // First and last terms

For i from 1 to n-1:
    x = a + i * h
    If i % 2 == 0:
        sum = sum + 2 * f(x) // Even index
    Else:
        sum = sum + 4 * f(x) // Odd index

Result = sum * (h / 3)
Return Result
```
**Further Study**
* [Simpson’s 1/3 Rule - GeeksforGeeks (Implementation)](https://www.geeksforgeeks.org/program-simpsons-13-rule/)
* [Simpson's Rule Derivation - Wolfram MathWorld](https://mathworld.wolfram.com/SimpsonsRule.html)

### 2. Simpson’s 3/8 Rule

**Theory: Cubic Approximation**
While the 1/3 rule uses parabolas (3 points), Simpson’s 3/8 Rule fits a **third-order polynomial (cubic curve)** through every four points. This generally provides slightly better accuracy for functions that are smoother.



The constraint for this method is that the number of intervals $n$ must be a **multiple of 3**.

$$I \approx \frac{3h}{8} \left[ (y_0 + y_n) + 3(y_1 + y_2 + y_4 + y_5 + \dots) + 2(y_3 + y_6 + \dots) \right]$$

**Algorithm**
1.  **Verify Intervals:** Ensure $n$ is a multiple of 3.
2.  **Calculate Step Size:** $h = (b - a) / n$.
3.  **Sum Extremes:** Add $y_0$ and $y_n$.
4.  **Sum Multiples of 3:** Multiply terms at indices divisible by 3 ($y_3, y_6...$) by 2.
5.  **Sum Others:** Multiply all remaining terms ($y_1, y_2, y_4, y_5...$) by 3.
6.  **Calculate Total:** Sum components and multiply by $3h/8$.

**Pseudocode**
```text
Input: Function f(x), lower_limit a, upper_limit b, intervals n
Output: Integral value

If n % 3 != 0:
    Print "Error: n must be divisible by 3"
    Return

h = (b - a) / n
sum = f(a) + f(b)

For i from 1 to n-1:
    x = a + i * h
    If i % 3 == 0:
        sum = sum + 2 * f(x) // Multiple of 3
    Else:
        sum = sum + 3 * f(x) // Rest of the terms

Result = sum * (3 * h / 8)
Return Result
```
**Further Study**
* [Simpson’s 3/8 Rule - GeeksforGeeks](https://www.geeksforgeeks.org/simpsons-38-rule-python/)
* [Numerical Integration Rules - Swarthmore College](https://lpsa.swarthmore.edu/NumInt/NumIntMain.html)
