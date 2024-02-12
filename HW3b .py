import math
# This code was made with help from ChatGPT
# This code was made with help from Jim Smay

# Simpson's rule for numerical integration
def simpsons_rule(f, a, b, n):
    """
    This function implements Simpson's rule for numerical integration.
    Args:
        f (function): The function to integrate.
        a (float): The lower limit of integration.
        b (float): The upper limit of integration.
        n (int): The number of intervals to use for the approximation.
    Returns:
        float: The approximate integral of f from a to b using Simpson's rule.
    """
    h = (b - a) / n
    s = f(a) + f(b)

    # Calculate the sum in the Simpson's rule formula
    for i in range(1, n, 2):
        s += 4 * f(a + i * h)
    for i in range(2, n-1, 2):
        s += 2 * f(a + i * h)

    return s * h / 3

# F-distribution function
def FZ(df, z):
    """
    This function calculates the value of the F-distribution function at z.
    Args:
        df (int): The degrees of freedom.
        z (float): The point at which to evaluate the function.
    Returns:
        float: The value of the F-distribution function at z.
    """
    # Define the probability density function of the F-distribution
    def pdf(t):
        return (math.gamma((df+1)/2) / (math.sqrt(df*math.pi) * math.gamma(df/2))) * (1 + t**2/df)**(-(df+1)/2)

    # Use Simpson's rule to calculate the integral of the pdf from -1000 to z
    return simpsons_rule(pdf, -1000, z, 10000)

def main():
    """
    This is the main function that prompts the user for the degrees of freedom and the upper limit of integration,
    then calculates and prints the value of the F-distribution function at the given point.
    """
    getOut = False
    while (getOut is False):
        m = int(input("Degrees of freedom (integer): "))
        u = float(input("Upper integration limit (float): "))
        result = FZ(m, u)
        print("F({:0.3f})={:0.3f}".format(u,result))
        getOut= input("Go again (Y/N)?") == "N"

if __name__ == "__main__":
    main()