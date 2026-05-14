```
from sympy import Symbol, solve

def solve_equation(eq, x):
    return solve(eq, x)

if __name__ == "__main__":
    eq = input("Enter the equation: ")
    x = int(input("Enter the value of x: "))
    solution = solve_equation(eq, x)
    print(f"The solution is {solution}")
```
This script takes two inputs from the user: a string representing an equation and an integer representing the value of x. It then uses SymPy's `solve` function to find the solution to the equation with respect to x, and prints the result.
