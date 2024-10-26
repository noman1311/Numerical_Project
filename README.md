## Abdullah Al Noman 2107105
## Avijit Ankon 2107083
## Abdullah Sheikh 2107035
# Numerical Methods Console Application

Our console application implements various numerical methods for solving linear equations, non-linear equations, differential equations, and matrix inversion.

## Table of Contents
- [Assignment Description](#assignment-description)
- [Application Structure](#application-structure)
  - [1. Solution of Linear Equations](#1-solution-of-linear-equations)
  - [2. Solution of Non-linear Equations](#2-solution-of-non-linear-equations)
  - [3. Solution of Differential Equations](#3-solution-of-differential-equations)
  - [4. Matrix Inversion](#4-matrix-inversion)
- [Technologies Used](#technologies-used)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)

---

## Assignment Description

This application is designed to solve various types of mathematical problems using numerical methods. It covers solutions for linear and non-linear equations, differential equations, and matrix operations eg. elementary matrix operations, matrix inversion. This project was completed in [Cpp].

## Application Structure

The application consists of solutions for linear equations, non-linear equations, differential equations, and matrix inversion, organized into 4 different header files. Specific functions were called inside one of the four header file to solve that particular problem. The structure is as follows-->
    1. What type of equation --> 2. Control goes to that header file. --> 3. Which method user want to use for solving that problem.

### 1. Solution of Linear Equations

Inside **linear_equation.h** header file.
Each method below provides a different approach to solving a system of linear equations:

#### a. Jacobi Iterative Method
The **Jacobi iterative method** is an algorithm used to solve a system of linear equations iteratively, particularly when the coefficient matrix is diagonally dominant. It’s useful in scenarios where direct methods (like Gaussian elimination) might be computationally expensive or infeasible for large systems.

### How It Works:
1. **Initialization**: Start with an initial guess for each variable, usually \( x = [0, 0, \dots, 0] \).
2. **Iteration**:
   - For each equation in the system, isolate the variable on the left side.
   - Calculate the new value of each variable based on the values from the previous iteration.
   - The formula for updating a variable \( x_i \) at iteration \( k+1 \) is:
     \[
     x_i^{(k+1)} = \frac{b_i - \sum_{j \neq i} a_{ij} x_j^{(k)}}{a_{ii}}
     \]
     where \( b_i \) is the constant term, \( a_{ij} \) are coefficients, and \( a_{ii} \) is the diagonal coefficient of \( x_i \).
3. **Convergence**:
   - Repeat the iterations until the difference between successive approximations is below a specified tolerance level, or until a maximum number of iterations is reached.
4. **Input**:
    - Matrix of coefficient, A.
    - Matrix of constants, B.
    - max tolerance and max iteration. 
5. **Output**:
    - A vector containing the results.
6. **Constraints**:
    - the matrix must be diagonally dominant. 
#### b. Gauss-Seidel Iterative Method
The **Gauss-Seidel method** is an iterative algorithm for solving a system of linear equations, similar to the Jacobi method but often with faster convergence. The key difference lies in how each variable is updated during the iteration.

### How It Works:
1. **Initialization**: Begin with an initial guess for each variable, typically \( x = [0, 0, \dots, 0] \).
2. **Iteration**:
   - For each variable \( x_i \), update its value immediately after calculating it within the same iteration.
   - This means that the updated value of \( x_i \) is used right away for calculating the subsequent variables within the same iteration.
   - The update formula for \( x_i \) at iteration \( k+1 \) is:
     \[
     x_i^{(k+1)} = \frac{b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)}}{a_{ii}}
     \]
     where \( b_i \) is the constant term, \( a_{ij} \) are coefficients, and \( a_{ii} \) is the diagonal coefficient for \( x_i \).

3. **Convergence**:
   - Continue iterating until the change between iterations is smaller than a specified tolerance or a maximum number of iterations is reached.
4. **Input**:
    - Matrix of coefficient, A.
    - Matrix of constants, B.
    - max tolerance and max iteration. 
5. **Output**:
    - A vector containing the results.
6. **Constraints**:
    - the matrix must be diagonally dominant.

#### c. Gauss Elimination
The **Gaussian elimination method** is a direct approach for solving a system of linear equations by transforming the coefficient matrix into an upper triangular form. This simplifies the system, allowing each variable to be solved sequentially through back-substitution.

### How It Works:
1. **Forward Elimination**:
   - The algorithm transforms the matrix into an upper triangular form by using row operations to eliminate the variables below the main diagonal.
   - Starting with the first row, each subsequent row is modified to make the entries below the pivot (diagonal element) zero by adding or subtracting multiples of the pivot row.

2. **Back Substitution**:
   - Once the matrix is in upper triangular form, start with the last equation and solve for the variable in that row.
   - Substitute this value back into the previous rows to find the remaining unknowns, moving up through the matrix.
3. **Input**:
•	Augmented Matrix (A): Coefficient matrix, size of n * (n+1).
4. **Output**:
The method outputs the solution vector (x) or the reduced row echelon form (RREF) of the augmented matrix.
5. **Corner Cases**:
•	If the coefficient matrix is singular (determinant is zero), the system may have no unique solution or infinitely many solutions.
•	If the matrix is poorly conditioned, numerical precision errors may lead to inaccurate results.

#### d. Gauss-Jordan Elimination
Gauss-Jordan elimination is an extension of Gauss Elimination that reduces the coefficient matrix to a diagonal form. This allows direct computation of each variable without requiring back-substitution, resulting in a straightforward solution for each unknown.
1. **Input**:
	Augmented Matrix (A): Coefficient matrix, size of n * (n+1).
2. **Output**:
The method outputs the solution vector (x) or the reduced row echelon form (RREF) of the augmented matrix.
3. **Corner Cases**:
•	If the coefficient matrix is singular (determinant is zero), the system may have no unique solution or infinitely many solutions.
•	If the matrix is poorly conditioned, numerical precision errors may lead to inaccurate results.

#### e. LU Factorization
LU Factorization decomposes the matrix of coefficients into a product of a lower triangular matrix \( L \) and an upper triangular matrix \( U \). Once the matrix is factorized, it simplifies solving the linear system by breaking it into two simpler systems that can be solved sequentially, especially useful when solving multiple systems with the same coefficient matrix.
1. **Input**:
•	Augmented Matrix (A): Coefficient matrix, size of n * (n+1).
2. **Output**:
•	The solution vector (x).
•	Decomposition of (A) into a lower triangular matrix (L) and an upper triangular matrix (U).
3. **Corner Cases**:
•	If matrix (A) is singular or ill-conditioned, LU factorization may fail.

**Note:** All methods are implemented for solving systems with a minimum of 5 equations.

### 2. Solution of Non-linear Equations

Under **non_linear_header.hpp** header file.
The following methods find approximate or exact root for non-linear equations \( f(x) = 0 \):

#### a. Bisection Method
The **Bisection method** is a numerical technique for finding roots of a continuous function \( f(x) = 0 \) within a given interval. It’s a bracketing method that works by repeatedly dividing an interval in half and selecting subintervals where a root must exist.

### How It Works:
1. **Initial Interval Selection**:
   - Start with an interval \([a, b]\) where the function values at \( a \) and \( b \) have opposite signs, indicating a root lies between them (by the Intermediate Value Theorem).
   
2. **Interval Halving**:
   - Calculate the midpoint \( c = \frac{a + b}{2} \).
   - Evaluate \( f(c) \):
     - If \( f(c) = 0 \), \( c \) is the root.
     - If \( f(a) \) and \( f(c) \) have opposite signs, set \( b = c \); otherwise, set \( a = c \).
   - Repeat this halving process until the interval is sufficiently small (defined by a tolerance).

3. **Convergence**:
   - The method converges slowly but reliably, with each iteration halving the interval size.
4. **Input**:
•	Vector of coefficient, tolerance, brackets
5. **Output**:
•	The solution variable (x).
6. **Constraints**:
    - f(a) * f(b) must be less than 0.  
    - **Guaranteed Convergence**: The method is reliable for continuous functions with opposite signs at the endpoints.
    - **Simple but Slow**: While very stable, the Bisection method can be slower than other methods, such as the Newton-Raphson method, particularly for high-precision requirements.

#### b. False Position Method
The False Position (Regula Falsi) method is a bracketing technique that improves upon the Bisection Method by calculating an approximation to the root using a weighted average between two points. This method can converge faster than the Bisection Method under certain conditions.

### How It Works:
1. **Initial Interval Selection**:
   - Start with an interval \([a, b]\) such that \( f(a) \) and \( f(b) \) have opposite signs, indicating the presence of a root.

2. **Secant Line Calculation**:
   - Calculate the point where the secant line (line connecting \( (a, f(a)) \) and \( (b, f(b)) \)) intersects the x-axis:
     \[
     c = b - \frac{f(b)(b - a)}{f(b) - f(a)}
     \]
   - Evaluate \( f(c) \):
     - If \( f(c) = 0 \), then \( c \) is the root.
     - If \( f(c) \) has the opposite sign of \( f(a) \), update \( b = c \); otherwise, update \( a = c \).

3. **Iteration**:
   - Repeat the process until \( f(c) \) is close enough to zero (within a defined tolerance).
4. **Input**:
•	Vector of coefficient, tolerance, brackets
5. **Output**:
•	The solution variable (x).
6. **Constraints**:
    - f(a) * f(b) must be less than 0.
    - **Faster Convergence**: Often converges more quickly than the Bisection method due to the secant approximation.
    - **Retains Bracketing**: Ensures that the root remains within the interval, providing reliability similar to the Bisection method. However, it may converge slowly if one endpoint is close to the actual root and does not change significantly during iterations.
#### c. Secant Method
The **Secant method** is an iterative technique for finding the root of a function \( f(x) = 0 \). Unlike bracketing methods (such as the Bisection or False Position methods), the Secant method uses two initial approximations and does not require the function values to have opposite signs. This approach is faster than the bracketing methods for functions that are well-behaved near the root.

### How It Works:
1. **Initial Guesses**:
   - Start with two initial approximations, \( x_0 \) and \( x_1 \), close to the suspected root of the function.

2. **Iteration Using Secant Line**:
   - Use the secant line between the points \( (x_0, f(x_0)) \) and \( (x_1, f(x_1)) \) to estimate the next point \( x_2 \) where \( f(x) \) is closer to zero:
     \[
     x_{n+1} = x_n - \frac{f(x_n)(x_n - x_{n-1})}{f(x_n) - f(x_{n-1})}
     \]
   - Replace \( x_{n-1} \) with \( x_n \), and \( x_n \) with \( x_{n+1} \), and repeat.

3. **Convergence**:
   - Continue the iterations until the difference between successive values of \( x \) is within a predefined tolerance, indicating that \( x_{n+1} \) is close enough to the root.
4. **Input**:
•	Vector of coefficient, tolerance, two assumptions
5. **Output**:
•	The solution variable (x).
#### d. Newton-Raphson Method
The **Newton-Raphson method** is a powerful and widely used iterative technique for finding the root of a function \( f(x) = 0 \). It is based on using the tangent line at a current approximation to estimate the next point closer to the root. The method converges rapidly for functions that are differentiable near the root, making it faster than most other methods.

### How It Works:
1. **Initial Guess**:
   - Start with an initial approximation \( x_0 \) close to the suspected root of \( f(x) \).

2. **Iteration Using Tangent Line**:
   - At each iteration, calculate the next approximation \( x_{n+1} \) using the formula:
     \[
     x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
     \]
   - Here, \( f'(x_n) \) is the derivative of \( f(x) \) at \( x_n \).

3. **Convergence**:
   - Repeat this process until \( |x_{n+1} - x_n| \) is smaller than a specified tolerance, indicating that \( x_{n+1} \) is sufficiently close to the root.
4. **Input**:
•	Vector of coefficient, tolerance, one assumption
5. **Output**:
•	The solution variable (x).
### 3. Solution of Differential Equations
  Under *ODE.h* header file
#### Runge-Kutta Method
The Runge-Kutta method (specifically the 4th-order version, often referred to as RK4) is a numerical technique for solving ordinary differential equations. It provides high accuracy by computing intermediate slopes within each interval, using a weighted average of these slopes to predict the solution at the next point. RK4 is widely used due to its balance of accuracy and computational efficiency.
**Input**:
•	Function (f): A function representing the first-order ODE of the form dy/dx = c1*x^(p1) + c2*y^(p2), where c1,c2,p1,p2 will be given by user.
 Initial condition (y0): The initial value of y.
•	Initial point (x0): The starting value of x.
•	Step size (h): The interval for computing values.
•	End point (x_end): The value of x at which the solution should be approximated.
**Output**:
•	Approximate values of y at different x values until reaching the endpoint.
**Corner Cases**:
•	Step size (h) that is too large may result in poor accuracy, while a very small step size increases computation time.
•	The function should be well-behaved over the interval for accurate results; discontinuities may cause issues.

### 4. Matrix Inversion
Under **matrix_inversion.h** header file

Matrix inversion is an operation that finds a matrix \( A^{-1} \) such that \( A \cdot A^{-1} = I \), where \( I \) is the identity matrix. This method is useful in solving matrix equations and systems of linear equations in general form. The implementation here involves transforming the matrix into an identity matrix through row operations.
**Input**:
•	Matrix (A): The square matrix to be inverted (n<=5).
**Output**:
•	The inverse of matrix (A), if it exists.
**Corner Cases**:
•	If the matrix (A) is singular (determinant is zero), it cannot be inverted.
•	For large matrices, numerical instability can lead to inaccurate results, especially when the matrix is poorly conditioned.

## Technologies Used

- [Cpp]
- GitHub for version control and collaboration

## Getting Started

### Prerequisites

- Install [Cpp]
- GitHub account

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/noman1311/Numerical_Project
<div align="center">
<table>

<tr>

  

<td  align="center">

<a  href="https://github.com/noman1311">

<img  src="https://scontent.fdac146-1.fna.fbcdn.net/v/t39.30808-6/428674404_2055657981488285_2769521409069830869_n.jpg?_nc_cat=107&ccb=1-7&_nc_sid=a5f93a&_nc_eui2=AeGE2ulnusk3pBDRCDuOA5OcVqe4VSYJ3CVWp7hVJgncJRnCkIUrZai2tO1zm4YLwjrOMhSKYQikY97Z1oOkDBbx&_nc_ohc=E3F52rXeInIQ7kNvgEuSE0_&_nc_pt=1&_nc_zt=23&_nc_ht=scontent.fdac146-1.fna&_nc_gid=Aaj-bWYp-BPtvuGtf0GuFEn&oh=00_AYB_D6NRBTENg0I4t_3sxQXZ0tJaLmMWWFxFSDtu_xBplQ&oe=6722DE81"  width="250px;"  alt="Abdullah Al Noman Profile Picture"/><br>

<sub>

<b>Abdullah Al Noman</b>

</sub>

</a>

</td>

  

<td  align="center">

<a  href="https://github.com/Avijit-35">

<img  src="https://scontent.fdac146-1.fna.fbcdn.net/v/t39.30808-6/345598524_202081039340833_2987918067526877551_n.jpg?_nc_cat=106&ccb=1-7&_nc_sid=6ee11a&_nc_eui2=AeGELphz6-ZDAIQJZmq9q_MZgh1ItO8HxTCCHUi07wfFMChAm6tWGnsIFJJN9ldudhV7mFrHsBEsUBPbICwal9Ma&_nc_ohc=Dnvp9V7p868Q7kNvgHhd7jE&_nc_pt=1&_nc_zt=23&_nc_ht=scontent.fdac146-1.fna&_nc_gid=A_aMLG0tPQi7SMP7Flq0qhl&oh=00_AYADV4bvL8jToHUZcASclybdfeFLNWi-Z-VlR_CsBAW-qg&oe=6722B36C"  width="250px;"  alt="Avijit Profile Picture"/><br>

<sub>

<b>Avijit Ankon</b>

</sub>

</a>

</td>

  

<td  align="center">

<a  href="https://github.com/crackhead6474">

<img  src="https://scontent.fdac146-1.fna.fbcdn.net/v/t39.30808-6/459418524_122105613728514030_4015433283834980416_n.jpg?_nc_cat=102&ccb=1-7&_nc_sid=6ee11a&_nc_eui2=AeGyfADJy6TnwT83np1Q0eDAC4CYEfo3jVILgJgR-jeNUrYjmHwMolOHjOoQ5uuXQt8CO7_p2oNJJLyV9xehU3VJ&_nc_ohc=M3I7BhC0bAYQ7kNvgHo6L6C&_nc_pt=1&_nc_zt=23&_nc_ht=scontent.fdac146-1.fna&_nc_gid=AEakkmO-hAd_zzG9i1ySKAZ&oh=00_AYDg9COfHprkqoG0N5N0UKSsupuRSPPeXau25uQErab9JA&oe=6722BE24"  width="250px;"  alt="Abdullah Sheikh Profile Picture"/><br>

<sub>

<b>Abdullah Sheikh</b>

</sub>

</a>

</td>
