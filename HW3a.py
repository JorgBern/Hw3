import DoolittleMethod as dm
import Gauss_Seidel as GS
import random
# ChatGPT was used to come up with thsi code
# Jim Smay helped with coming up this code by providing Doolittle and Gauss_seidel

# Cholesky decomposition method for solving linear systems
def Cholesky(Aaug):
    """
    This function performs the Cholesky decomposition on a given matrix.
    It first separates the augmented matrix into A and b, then constructs the lower triangular matrix L.
    It then solves the system L*y = b for y and L^T*x = y for x.
    Args:
        Aaug (list): The augmented matrix [A|b]
    Returns:
        x (list): The solution vector
        L (list): The lower triangular matrix
        Ltrans (list): The transpose of the lower triangular matrix
    """
    # Separate the augmented matrix into A and b
    A, b = GS.separateAugmented(Aaug)
    n= len(A)
    # Initialize the lower triangular matrix L
    L = [[0] * n for _ in range(n)]

    # Construct the lower triangular matrix L
    for i in range(n):
        for j in range(i + 1):
            if i == j:
                temp_sum = sum(L[i][k] ** 2 for k in range(j))
                L[i][j] = (A[i][j] - temp_sum) ** 0.5
            else:
                temp_sum = sum(L[i][k] * L[j][k] for k in range(j))
                L[i][j] = (A[i][j] - temp_sum) / L[j][j]

    # Compute the transpose of L
    Ltrans = transpose(L)

    # Solve the system L*y = b for y using forward substitution
    y = [0 for i in range(n)]
    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))
        y[i] = y[i] / L[i][i]

    # Solve the system L^T*x = y for x using back substitution
    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        x[i] = y[i] - sum(Ltrans[i][j] * x[j] for j in range(i+1, n))
        x[i] = x[i] / Ltrans[i][i]

    return x, L, Ltrans

# Check if a matrix is symmetric positive definite
def SymPosDef(A):
    """
    This function checks if a given matrix is symmetric positive definite.
    It first computes the transpose of the matrix, then checks if the matrix is symmetric.
    It then generates a random vector x and checks if x_transpose * A * x is positive.
    Args:
        A (list): The matrix to check
    Returns:
        bool: True if the matrix is symmetric positive definite, False otherwise
    """
    # Compute the transpose of A
    A_transpose = transpose(A)

    # Check if A is symmetric
    for i in range(len(A)):
        for j in range(len(A[0])):
            if A[i][j] != A_transpose[i][j]:
                return False

    # Generate a random vector x
    n = len(A)
    x = [random.uniform(-1, 1) for _ in range(n)]

    # Compute x_transpose * A * x
    x_trans_A = sum(x[i] * sum(A[i][j] * x[j] for j in range(n)) for i in range(n))

    # Check if x_transpose * A * x is positive
    if x_trans_A > 0:
        return True
    else:
        return False

# Transpose a matrix
def transpose(A):
    """
    This function computes the transpose of a given matrix.
    Args:
        A (list): The matrix to transpose
    Returns:
        list: The transpose of the matrix
    """
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]

def main():
    """
    This is the main function that defines the given matrices and solves the linear systems using either the Cholesky or Doolittle method.
    """
    # Define the given matrices
    A1_aug = [[1, -1, 3, 2, 15], [-1, 5, -5, -2, -35], [3, -5, 19, 3, 94], [2, -2, 3, 21, 1]]
    A2_aug = [[4, 2, 4, 0, 20], [2, 2, 3, 2, 36], [4, 3, 6, 3, 60], [0, 2, 3, 9, 122]]

    # Solve the linear systems for each matrix
    for idx, Aaug in enumerate([A1_aug, A2_aug], start=1):
        # Separate the augmented matrix into A and b
        A, b = GS.separateAugmented(Aaug)

        # Check if A is symmetric positive definite
        if SymPosDef(A):
            # If A is symmetric positive definite, use the Cholesky method
            x, _, _ = Cholesky(Aaug)
            method = "Cholesky"
        else:
            # If A is not symmetric positive definite, use the Doolittle method
            x = dm.Doolittle(Aaug)
            method = "Doolittle"

        # Print the solution vector and the method used
        print(f"Solution vector for Matrix {idx} using {method} method:", x)
        print()

if __name__ == "__main__":
    main()
