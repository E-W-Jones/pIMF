# Create some methods for dealing with interpolated data/linear piecewise functions
def interpolate(t, x, y, N=None, extrapolate=False):
    """For data defined (xi, yi), calculate y(t). Optionally provide the length, N, and extrapolation."""
    if N is None:
        N = len(x)

    if x[1] < x[0]:
        raise ValueError("It looks like your x array is not sorted.")

    for i in range(N-1):
        if x[i] <= t <= x[i+1]:
            return y[i] + (t - x[i]) * (y[i+1] - y[i]) / (x[i+1] - x[i])

    if extrapolate is True:
        if t < x[0]:
            return y[0] + (t - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
        elif t > x[-1]:
            return y[-2] + (t - x[-2]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    else:
        raise ValueError(f"{t} lies outwith the bounds of ({x[0]}, {x[-1]})")

def find_bracketing_indices(t, x, N=None):
    """Find the indices i, j st xi <= t < xj. If t is outside of min(x), max(x), return edge index twice."""
    if N is None:
        N = len(x)
    # Just go left to right increasing.
    if t <= x[0]:
        return 0, 0
    for i in range(N-1):
        if x[i] <= t < x[i+1]:
            return i, i+1
    else:
        return i+1, i+1

def integrate(a, b, x, y, N=None, debug=False):
    """Integrate the linear piecewise function defined bt xi, yi from x=a to x=b."""
    #TODO: What's up with those boundary conditions huh?
    if N is None:
        N = len(x)
    ai = bi = None
    for i in range(N-1):
        if x[i] <= a < x[i+1]:
            ai = i
        if x[i] <= b < x[i+1]:
            bi = i
    if a == x[N-1]:
        ai = N-2
    if b == x[N-1]:
        bi = N-2

    if debug:
        print(f"{ai = }, {bi = }")

    if ai is None:
        raise ValueError("ai is None. Might currently be where I don't have boundary conditions.")
    if bi is None:
        raise ValueError("bi is None. Might currently be where I don't have boundary conditions.")

    if ai == bi:
        i = ai
        m = (y[i+1] - y[i]) / (x[i+1] - x[i])
        return (y[i] - m*x[i] + m*(b+a)/2) * (b - a)

    # Otherwise we need to calculate from a -> ai+1, ai+1 -> bi, and bi -> b seperately
    integral = 0
    # a-> ai+1
    m = (y[ai+1] - y[ai]) / (x[ai+1] - x[ai])
    integral += (y[ai] - m*x[ai] + m*(x[ai+1]+a)/2) * (x[ai+1] - a)
    # bi -> b
    m = (y[bi+1] - y[bi]) / (x[bi+1] - x[bi])
    integral += (y[bi] - m*x[bi] + m*(b + x[bi])/2) * (b - x[bi])
    # ai+1 -> bi
    for i in range(ai+1, bi):
        print(f"{i = }")
        integral += 0.5 * (x[i+1] - x[i]) * (y[i+1] + y[i])
    return integral

def integrate_weighted(a, b, x, y, weight_integral, weight_integral_product, N=None, debug=False):
    """
    Calculate the integral if a piecewise linear function defined by (x, y) from a to b, where the integrand is weighted by some function f
    int^b_a f(x) y(x) dx

    Needs weight_integral, which should be a function that takes a and b as paramaters, returning int^b_a f(x)dx,
    and weight_integral_product, which is similar, but the integrand is the product of x and f(x) is int^b_a xf(x)dx
    """
    if N is None:
        N = len(x)
    ai = bi = None
    for i in range(N-1):
        if x[i] <= a < x[i+1]:
            ai = i
        if x[i] <= b < x[i+1]:
            bi = i
    if a == x[N-1]:
        ai = N-2
    if b == x[N-1]:
        bi = N-2

    if debug:
        print(f"{ai = }, {bi = }")


    if ai is None:
        print(a, b, x, y)
        raise ValueError("ai is None. Might currently be where I don't have boundary conditions.")
    if bi is None:
        raise ValueError("bi is None. Might currently be where I don't have boundary conditions.")

    if ai == bi:
        i = ai
        m = (y[i+1] - y[i]) / (x[i+1] - x[i])
        return (y[i] - m*x[i]) * weight_integral(a, b) + m * weight_integral_product(a, b)

    # Otherwise we need to calculate from a -> ai+1, ai+1 -> bi, and bi -> b seperately
    integral = 0
    # a-> ai+1
    m = (y[ai+1] - y[ai]) / (x[ai+1] - x[ai])
    integral += (y[ai] - m*x[ai]) * weight_integral(a, x[ai+1]) + m * weight_integral_product(a, x[ai+1])
    # bi -> b
    m = (y[bi+1] - y[bi]) / (x[bi+1] - x[bi])
    integral += (y[bi] - m*x[bi]) * weight_integral(x[bi], b) + m * weight_integral_product(x[bi], b)
    # ai+1 -> bi
    for i in range(ai+1, bi):
        if debug:
            print(f"{i = }")
        m = (y[i+1] - y[i]) / (x[i+1] - x[i])
        integral += (y[i] - m*x[i]) * weight_integral(x[i], x[i+1]) + m * weight_integral_product(x[i], x[i+1])
    return integral

def test_integrate(a=None, b=None, x=None, y=None, debug=True):
    import numpy as np
    from scipy.integrate import quad
    import matplotlib.pyplot as plt

    if x is None:
        x = np.array([0, 1, 2, 3])
    if y is None:
        y = x**2
    if a is None:
        a = min(x)
    if b is None:
        b = max(x)


    print(integrate(a, b, x, y, debug=debug))
    print(quad(lambda x_: np.interp(x_, x, y), a, b))

def test_integrate_product(a=None, b=None, x=None, y=None, weight=None, weight_integral=None, weight_integral_product=None, debug=True):
    import numpy as np
    from scipy.integrate import quad
    import matplotlib.pyplot as plt

    if x is None:
        x = np.array([0, 1, 2, 3])
    if y is None:
        y = x**2
    if a is None:
        a = min(x)
    if b is None:
        b = max(x)
    if weight is None:
        def weight(x_):
            return 1
    if weight_integral is None:
        def weight_integral(a_, b_):
            return quad(weight, a_, b_)[0]
    if weight_integral_product is None:
        def weight_integral_product(a_, b_):
            return quad(lambda x: x * weight(x), a_, b_)[0]

    print(integrate_weighted(a, b, x, y, weight_integral, weight_integral_product, debug=debug))
    print(quad(lambda x_: weight(x_) * np.interp(x_, x, y), a, b))


def main():
    import numpy as np
    print("====== Testing interpolate ======")
    x = np.array([0, 1, 2, 3, 4])
    y = x**2
    for t, yexpected in zip([-1, 0, 1.5, 4, 5],
                            [-1, 0, 3*1.5 - 2, 16, 7*5 - 12]
                            ):  # Expected for extrapolating linearly at edges
        print(f"Me: {interpolate(t, x, y, extrapolate=True)} expected: {yexpected} diff: {interpolate(t, x, y, extrapolate=True) - yexpected}")

    try:
        interpolate(-1, x, y)
    except ValueError as e:
        print(f"Successfully raised error with message: {e.args}")

    print("====== Testing integrate ======")
    test_integrate(0, 0, x=np.array([0, 1, 2, 3]))  # ai = 0, bi = 0
    test_integrate(1, 1.5, x=np.array([0, 1, 2, 3]), y=np.array([0, 1, 2, 3]))  # ai = 0, bi = 0
    test_integrate(1, 1.5, x=np.array([0, 1, 2, 3]), y=np.array([0, 1, 2, 3])**2)  # ai = 0, bi = 0
    test_integrate(3, 3, x=np.array([0, 1, 2, 3]))  # ai = 3, bi = 3
    test_integrate(0.5, 2.5, x=np.array([0, 1, 2, 3]))  # ai = 0, bi = 2

    print("====== Testing integrate product ======")
    test_integrate_product(0, 0)
    test_integrate_product(3, 3)
    test_integrate_product(0.5, 2.5)
    test_integrate_product(0, 3)

    def weight(x):
        return x

    print("=== Different weight ===")
    test_integrate_product(0, 0, weight=weight)
    test_integrate_product(3, 3, weight=weight)
    test_integrate_product(0.5, 2.5, weight=weight)
    test_integrate_product(0, 3, weight=weight)

    print("=== With defined weight integrals ===")

    def weight_integral(a, b):
        return (b*b - a*a) / 2

    def weight_integral_product(a, b):
        return (b*b*b - a*a*a) / 3

    test_integrate_product(0, 0, weight=weight, weight_integral=weight_integral, weight_integral_product=weight_integral_product)
    test_integrate_product(3, 3, weight=weight, weight_integral=weight_integral, weight_integral_product=weight_integral_product)
    test_integrate_product(0.5, 2.5, weight=weight, weight_integral=weight_integral, weight_integral_product=weight_integral_product)
    test_integrate_product(0, 3, weight=weight, weight_integral=weight_integral, weight_integral_product=weight_integral_product)

if __name__ == "__main__":
    main()