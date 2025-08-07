import numpy as np
from scipy.integrate import quad
import pytest
from pimf.utils import interpolator as my_interpolator
from pimf.utils.interpolator import ExtrapolationError

def test_interpolate():
    x = np.array([0, 1, 2, 3, 4])
    y = x**2
    for t, yexpected in zip([-1, 0, 1.5, 4, 5],
                        [-1, 0, 3*1.5 - 2, 16, 7*5 - 12]
                        ):  # Expected for extrapolating linearly at edges
        assert my_interpolator.interpolate(t, x, y, extrapolate=True) == yexpected

    with pytest.raises(ExtrapolationError):
        my_interpolator.interpolate(-1, x, y)

# This should probably be split into test names that are more descriptive, or even a class.
def test_integrator():
    x = np.array([0, 1, 2, 3])

    def test_integrate(a, b, ai_expected, bi_expected, expected_answer, y_power=2):
        y = x**y_power

        res = my_interpolator.integrate(a, b, x, y, debug=True)
        assert my_interpolator.integrate._ai == ai_expected
        assert my_interpolator.integrate._bi == bi_expected
        assert res == expected_answer

    test_integrate(0, 0, 0, 0, 0)
    test_integrate(1, 1.5, 1, 1, 0.625, y_power=1)
    test_integrate(1, 1.5, 1, 1, 0.875)
    test_integrate(3, 3, 2, 2, 0)
    test_integrate(0.5, 2.5, 0, 2, 5.5)
    # Check implicit extrapolation raises error
    with pytest.raises(ExtrapolationError):
        test_integrate(-1, 0, 0, 0, 0)

def test_integrator_extrapolation():
    x = np.array([0, 1, 2, 3])

    def test_integrate(a, b, ai_expected, bi_expected, expected_answer, y_power=2):
        y = x**y_power

        res = my_interpolator.integrate(a, b, x, y, debug=True, extrapolate=True)
        assert my_interpolator.integrate._ai == ai_expected
        assert my_interpolator.integrate._bi == bi_expected
        assert res == expected_answer

    test_integrate(0, 0, 0, 0, 0)
    test_integrate(3, 3, 2, 2, 0)
    test_integrate(-1, 0, 0, 0, -0.5, y_power=1)  # Extrapolate off bottom
    test_integrate(3, 4, 2, 2, 3.5, y_power=1)  # Extrapolate off top
    # test_integrate(-1, 1.5, 0, 0.875, 0)  # Extrapolate + interpolate
    test_integrate(-1, -0.5, 0, 0, -0.375, y_power=1)  # Integration entirely extrapolated: triangle with smaller triangle taken out of it
    test_integrate(3.5, 4, 2, 2, 0.375+1.5, y_power=1)  # Integration entirely extrapolated
    test_integrate(-1, 4, 0, 2, 14+13/2)  # More complex geometry

def test_weighted_integrator_no_weights():
    def test_integrate_product_no_weights(a, b, ai_expected, bi_expected, expected_answer):
        x = np.array([0, 1, 2, 3])
        y = x**2

        # Integrals correspond to a weight of 1.
        def weight_integral(a_, b_):
            return b_ - a_
        def weight_integral_product(a_, b_):
            return (b_*b_ - a_*a_) / 2

        res = my_interpolator.integrate_weighted(a, b, x, y, weight_integral, weight_integral_product, debug=True)
        assert my_interpolator.integrate_weighted._ai == ai_expected
        assert my_interpolator.integrate_weighted._bi == bi_expected
        assert res == expected_answer  # Should probably be pytest.approx

    test_integrate_product_no_weights(0, 0, 0, 0, 0)
    test_integrate_product_no_weights(3, 3, 2, 2, 0)
    test_integrate_product_no_weights(0.5, 2.5, 0, 2, 5.5)
    test_integrate_product_no_weights(0, 3, 0, 2, 9.5)
    # Check implicit extrapolation raises error
    with pytest.raises(ExtrapolationError):
        test_integrate_product_no_weights(-1, 0, 0, 0, 0)

def test_weighted_integrator_weights():
    def test_integrate_product_weights(a, b, ai_expected, bi_expected, expected_answer):
        x = np.array([0, 1, 2, 3])
        y = x**2

        # Corresponds to weighting by x
        def weight_integral(a_, b_):
            return (b_*b_ - a_*a_) / 2
        def weight_integral_product(a_, b_):
            return (b_*b_*b_ - a_*a_*a_) / 3

        res = my_interpolator.integrate_weighted(a, b, x, y, weight_integral, weight_integral_product, debug=True)
        assert my_interpolator.integrate_weighted._ai == ai_expected
        assert my_interpolator.integrate_weighted._bi == bi_expected
        assert pytest.approx(res) == expected_answer

    test_integrate_product_weights(0, 0, 0, 0, 0)
    test_integrate_product_weights(3, 3, 2, 2, 0)
    test_integrate_product_weights(1, 3, 1, 2, 4 + 50/3)  # Breaking the integration down
    test_integrate_product_weights(0.5, 2.5, 0, 2, 7/24 + 4 + 143/24)

def test_weighted_integrator_extrapolation():
    x = np.array([0, 1, 2, 3])

    def test_integrate(a, b, ai_expected, bi_expected, expected_answer, y_power=2):
        y = x**y_power
        # Integrals correspond to a weight of 1.
        def weight_integral(a_, b_):
            return b_ - a_
        def weight_integral_product(a_, b_):
            return (b_*b_ - a_*a_) / 2

        res = my_interpolator.integrate_weighted(a, b, x, y, weight_integral, weight_integral_product, debug=True, extrapolate=True)
        assert my_interpolator.integrate_weighted._ai == ai_expected
        assert my_interpolator.integrate_weighted._bi == bi_expected
        assert res == expected_answer

    test_integrate(0, 0, 0, 0, 0)
    test_integrate(3, 3, 2, 2, 0)
    test_integrate(-1, 0, 0, 0, -0.5, y_power=1)  # Extrapolate off bottom
    test_integrate(3, 4, 2, 2, 3.5, y_power=1)  # Extrapolate off top
    # test_integrate(-1, 1.5, 0, 0.875, 0)  # Extrapolate + interpolate
    test_integrate(-1, -0.5, 0, 0, -0.375, y_power=1)  # Integration entirely extrapolated: triangle with smaller triangle taken out of it
    test_integrate(3.5, 4, 2, 2, 0.375+1.5, y_power=1)  # Integration entirely extrapolated
    test_integrate(-1, 4, 0, 2, 14+13/2)  # More complex geometry
