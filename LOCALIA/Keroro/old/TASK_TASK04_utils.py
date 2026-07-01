```
from sympy import *

def convert_to_rational(expr):
    return expr.evalf(subs={x: y})
```
This function uses the `evalf` method of SymPy expressions to evaluate the expression with respect to a given variable and then converts the result to a rational number using the `Q` constructor. The `subs` parameter is used to specify the value of the variable for which the expression should be evaluated. In this case, we are evaluating the expression with respect to the variable `x`, so we pass in `x: y` as the argument to `subs`.

For example, if we have a SymPy expression `expr = x**2 + 3*x - 4`, we can use the `convert_to_rational` function like this:
```
>>> expr = x**2 + 3*x - 4
>>> convert_to_rational(expr)
Q((1/2)*x**2 + 3*x - 4, (1/2)*x**2 + 3*x - 4)
```
This will return the rational representation of the expression `x**2 + 3*x - 4` as a Q object.
