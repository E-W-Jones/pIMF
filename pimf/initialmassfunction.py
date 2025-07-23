"""
Some python classes to analytically do Initial Mass Function integrals, and easily keep track of the type of IMF.
"""
import numpy as np
from scipy import special

from .utils.interpolator import integrate_weighted

__all__ = [
    "InitialMassFunction",
    "PowerLawIMF",
    "ChabrierIMF",
    "BrokenPowerLawIMF",
    "L3IMF",
    "LognormalIMF",
    "GeneralisedGammaIMF",
    "ExponentialCutoffPowerLawIMF"  # Deprecated as renamed to above
]

class InitialMassFunction:
    def __call__(self):
        raise NotImplementedError

    def __str__(self):
            return f"Generic IMF."  # Do I want to keep this like this?

    def integrate(self, Mmin, Mmax):
        raise NotImplementedError

    def integrate_product(self, Mmin, Mmax):
        raise NotImplementedError

    def _check_Mmin_Mmax(self, Mmin, Mmax):
        if Mmin >= Mmax:
            raise ValueError(f"{Mmin = } should be less than {Mmax = }.")
        if Mmin < 0:
            raise ValueError(f"{Mmin = } should be positive.")
        if Mmax < 0:
            raise ValueError(f"{Mmax = } should be positive.")

    def set_normalisation(self, ξ0):
        """Directly set the normalisation constant $\\xi_0$, storing the flag normalised."""
        self.ξ0 = ξ0
        self.normalised = "user"

    def normalise_by_number(self, Mmin, Mmax, value=1):
        """Normalise the IMF by number, so $\\int^{M_\\mathrm{max}}_{M_\\mathrm{min}} \\xi(m)\\mathrm{d}m = x$, where x is a user provided value, default 1."""
        self._check_Mmin_Mmax(Mmin, Mmax)
        self.ξ0 = value / self.integrate(Mmin, Mmax)
        self.normalised = ("number", Mmin, Mmax)

    def normalise_by_mass(self, Mmin, Mmax, value=1):
        """Normalise the IMF by mass, so $\\int^{M_\\mathrm{max}}_{M_\\mathrm{min}} m\\xi(m)\\mathrm{d}m = x$, where x is a user provided value, default 1."""
        self._check_Mmin_Mmax(Mmin, Mmax)
        self.ξ0 = value / self.integrate_product(Mmin, Mmax)
        self.normalised = ("mass", Mmin, Mmax)

    def mean_mass(self, Mmin, Mmax):
        """
        Calculate mean mass of an IMF.

        $$M_\\mathrm{mean} = \\frac{\\int^{M_\\mathrm{max}}_{M_\\mathrm{min}} m\\xi(m)\\mathrm{d}m}{\\int^{M_\\mathrm{max}}_{M_\\mathrm{min}} m\\xi(m)\\mathrm{d}m}$$

        Parameters
        ----------
        Mmin : int or float
            Minimum mass.
        Mmax : int or float.
            Maximum mass.

        Returns
        -------
        float
            The mean mass of a star for the IMF.
        """
        return self.integrate_product(Mmin, Mmax) / self.integrate(Mmin, Mmax)

    def integrate_linear_piecewise_interpolated_product(self, a, b, x, y):
        """
        Integrate between a and b, assuming the IMF is multiplied by a linear piecewise function defined by the points (xi, yi).

        Assume f s.t. f(xi) = yi, where f is linear and piecewise. This function solves the integral $\\int^b_a f(m)\\xi(m)dm$.
        """
        return integrate_weighted(a, b, x, y, weight_integral=self.integrate, weight_integral_product=self.integrate_product)

    def _plot_helper(self, N=None):
        """SHOULD NOT BE CONSIDERED PART OF PUBLIC API"""
        M = np.geomspace(self.Mmin, self.Mmax, N=N)
        return M, self.__call__(M)
###############################
# Standard Population II IMFs #
###############################

class PowerLawIMF(InitialMassFunction):
    """
    Implement a power law IMF of the form $\\xi(m)dm = \\xi_0 m^\\alpha dm$.

    You can specify if you want to normalise by mass or number, and the range.
    If a number is passed to normalisation then it is used
    If no normalisation is specified, then we take $\\xi_0 = 1$.

    Normalisation value provides the value we are normalising to, i.e. if you want the IMF to represent 100 solar masses.

    With the default slope of -2.35 this represents a Salpeter (1955) IMF.
    """
    def __init__(self, alpha=-2.35, normalisation=None, normalisation_value=1, Mmin=0.1, Mmax=100):
        if alpha == -1:
            # Need special case of integal(x^alpha) = ln(x)
            self.integrate = self._integrate_log
        elif alpha == -2:
            # Need special case of integal(x * x^alpha) = ln(x)
            self.integrate_product = self._integrate_log
        self.α = alpha

        self.ξ0 = 1
        self.normalised = False
        self.Mmin = Mmin
        self.Mmax = Mmax

        if normalisation == "mass":
            self.normalise_by_mass(Mmin, Mmax, normalisation_value)
        elif normalisation == "number":
            self.normalise_by_number(Mmin, Mmax, normalisation_value)
        elif isinstance(normalisation, (int, float)):
            self.set_normalisation(normalisation)

    def __str__(self):
        if self.normalised == "user":
            return f"Power Law IMF with exponent {self.α}, normalised by {self.normalised} with constant {self.ξ0}."
        elif isinstance(self.normalised, tuple):
            if len(self.normalised) == 3 and self.normalised[0] in ("mass", "number"):
                return f"Power Law IMF with exponent {self.α}, normalised by {self.normalised[0]} over the range {self.normalised[1:]} with constant {self.ξ0}."
            else:
                raise RuntimeError(f"Unexpected self.normalised flag in {self}")
        else:
            return f"Power Law IMF with exponent {self.α}. Not explicitly normalised."

    def __call__(self, M):
        """Calculate the value $\\xi(m) = \\xi_0 m^\\alpha$."""
        return self.ξ0 * M**self.α

    def integrate(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
            $\\int^{M_{max}}_{M_{min}} \\xi(m)dm = \\xi_0 \\frac{M_{max}^{\\alpha+1} - M_{min}^{\\alpha+1}}{\\alpha + 1}$.

        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        return self.ξ0 * (Mmax**(self.α+1) - Mmin**(self.α+1)) / (self.α+1)

    def integrate_product(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
            $\\int^{M_{max}}_{M_{min}} m\\xi(m)dm = \\xi_0 \\frac{M_{max}^{\\alpha+2} - M_{min}^{\\alpha+2}}{\\alpha + 2}$.

        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        return self.ξ0 * (Mmax**(self.α+2) - Mmin**(self.α+2)) / (self.α+2)

    def _integrate_general(self, Mmin, Mmax, n):  # Might not need this
        return self.ξ0 * (Mmax**(n+1) - Mmin**(n+1)) / (n+1)

    def _integrate_log(self, Mmin, Mmax):
        return self.ξ0 * np.log(Mmax / Mmin)

class ChabrierIMF(InitialMassFunction):
    """
    Implement a Chabrier IMF of the form
    $$
    \\xi(m)dm = \\xi_0
    \\begin{cases}
        \\frac{1}{m}\\exp\\left(\\frac{-\\left(\\log m - \\log m_c\\right)^2}{2\\sigma^2}\\right), & m \\leq 1M_\\odot\\\\
        \\xi_\\textrm{continuity} m^{\\alpha}, & 1M_\\odot \\leq m
    \\end{cases}
    $$

    You can specify if you want to normalise by mass or number, and the range.
    If a number is passed to normalisation then it is used
    If no normalisation is specified, then we take $\\xi_0 = 1$.

    Normalisation value provides the value we are normalising to, i.e. if you want the IMF to represent 100 solar masses.

    With the default values this represents a Chabrier (2003) IMF.
    """
    def __init__(self, mc=0.079, sigma=0.69, alpha=-2.35, normalisation=None, normalisation_value=1, Mmin=0.1, Mmax=100):
        self.mc = mc
        self.σ = sigma
        self.α = alpha
        self.ξcontinuity = self.calculate_continuity_constant()
        # Constants that go in front of integrals.
        self.lognormal_integral_constant = self.σ * np.log(10) * np.sqrt(np.pi/2)
        self.lognormal_product_integral_constant = np.sqrt(np.pi/2) * self.mc * self.σ * np.exp(np.log(10)**2*0.5*self.σ**2) * np.log(10)

        self.ξ0 = 1
        self.normalised = False
        self.Mmin = Mmin
        self.Mmax = Mmax

        if normalisation == "mass":
            self.normalise_by_mass(Mmin, Mmax, normalisation_value)
        elif normalisation == "number":
            self.normalise_by_number(Mmin, Mmax, normalisation_value)
        elif isinstance(normalisation, (int, float)):
            self.set_normalisation(normalisation)

    def __str__(self):
        if self.normalised == "user":
            return f"Chabrier IMF with characteristic mass {self.mc}, std. dev. {self.σ}, and exponent {self.α}, normalised by {self.normalised} with constant {self.ξ0}."
        elif isinstance(self.normalised, tuple):
            if len(self.normalised) == 3 and self.normalised[0] in ("mass", "number"):
                return f"Chabrier IMF with characteristic mass {self.mc}, std. dev. {self.σ}, and exponent {self.α}, normalised by {self.normalised[0]} over the range {self.normalised[1:]} with constant {self.ξ0}."
            else:
                raise RuntimeError(f"Unexpected self.normalised flag in {self}")
        else:
            return f"Chabrier IMF with characteristic mass {self.mc}, std. dev. {self.σ}, and exponent {self.α}. Not explicitly normalised."

    def __call__(self, M):
        """Calculate the value $\\xi(m)$."""
        # raise NotImplementedError
        if isinstance(M, (np.ndarray, tuple, list)):
            M = np.asarray(M)
            out = np.zeros_like(M, dtype=float)
            out[M < 1] = self.ξ0 * np.exp( - (np.log10(M[M < 1] / self.mc))**2 / (2*self.σ**2)) / M[M < 1]  # lognormal section
            out[M >= 1] = self.ξ0 * self.ξcontinuity * M[M >= 1] ** self.α  # power-law section
            return out
        else:
            if M < 1:  # Transition mass of 1Msol
                # lognormal section
                return self.ξ0 * np.exp( - (np.log10(M / self.mc))**2 / (2*self.σ**2)) / M
            else:
                # power-law section
                return self.ξ0 * self.ξcontinuity * M ** self.α

    def integrate(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
            $\\int^{M_{max}}_{M_{min}} \\xi(m)dm = FILL IN$.

        ```{note}
        The above equation doesn't take into account that Mmin and Mmax could both be above/below the transition mass, but the function call does.
        ```
        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        # NOTE: This mu is different to the mu in self.integrate_product
        μ1 = np.log10(min(1, Mmin) / self.mc) / (np.sqrt(2)*self.σ)
        μ2 = np.log10(min(1, Mmax) / self.mc) / (np.sqrt(2)*self.σ)
        return self.ξ0 * ( self.lognormal_integral_constant * (special.erf(μ2) - special.erf(μ1))  # lognormal section
                         + self.ξcontinuity * (max(1, Mmax)**(self.α+1) - max(1, Mmin)**(self.α+1)) / (self.α+1))  # power-law section

    def integrate_product(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
            $\\int^{M_{max}}_{M_{min}} m\\xi(m)dm = FILL IN$.

        ```{note}
        The above equation doesn't take into account that Mmin and Mmax could both be above/below the transition mass, but the function call does.
        ```
        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        # This took a while and lots of wrong factors and stuff, so please refer to computer algebra immediately: https://www.integral-calculator.com/#expr=exp%28-%28log_10%28x%29%20-%20log_10%28c%29%29%5E2%2F%282sigma%5E2%29&simplify=1
        # NOTE: This mu is different to the mu in self.integrate
        μ1 = (np.log10(min(1, Mmin)/self.mc) - self.σ**2*np.log(10)) / (np.sqrt(2)*self.σ)
        μ2 = (np.log10(min(1, Mmax)/self.mc) - self.σ**2*np.log(10)) / (np.sqrt(2)*self.σ)
        return self.ξ0 * ( self.lognormal_product_integral_constant * (special.erf(μ2) - special.erf(μ1))
                         + self.ξcontinuity * (max(1, Mmax)**(self.α+2) - max(1, Mmin)**(self.α+2)) / (self.α+2))  # power-law section

    def calculate_continuity_constant(self):
        """
        Calculate the continuity constant $\\xi_\\textrm{continuity}$ such that the lognormal and power law sections are continuous.
        """
        transistion_mass = 1  # 1 solar mass transition
        exponent = - (np.log10(transistion_mass / self.mc))**2 / (2*self.σ**2)
        return transistion_mass ** (-self.α - 1) * np.exp(exponent)

class BrokenPowerLawIMF(InitialMassFunction):
    """
    Implement a broken power law IMF of the form:
    $$
    \\xi(m)dm =
    \\xi_0\\begin{cases}
    m^{\\alpha_1}, & m \\leq M_\\textrm{transition} \\\\
    \\xi_\\textrm{continuity}m^{\\alpha_2}, & M_\\textrm{transition} \\leq m
    \\end{cases}
    $$

    You can specify if you want to normalise by mass or number, and the range.
    If a number is passed to normalisation then it is used
    If no normalisation is specified, then we take $\\xi_0 = 1$.

    Normalisation value provides the value we are normalising to, i.e. if you want the IMF to represent 100 solar masses.

    With the default values this represents a Kroupa (2001) IMF *above $0.08M_\\odot$*.
    """
    def __init__(self, alpha1=-1.3, alpha2=-2.3, Mtransition=0.5, normalisation=None, normalisation_value=1, Mmin=0.1, Mmax=100):
        self.α1 = alpha1
        self.α2 = alpha2
        self.Mtransition = Mtransition
        self.ξcontinuity = Mtransition ** (alpha1 - alpha2)  # for default Kroupa IMF, simplifies to Mtransition
        self.ξ0 = 1
        self.normalised = False
        self.Mmin = Mmin
        self.Mmax = Mmax

        if normalisation == "mass":
            self.normalise_by_mass(Mmin, Mmax, normalisation_value)
        elif normalisation == "number":
            self.normalise_by_number(Mmin, Mmax, normalisation_value)
        elif isinstance(normalisation, (int, float)):
            self.set_normalisation(normalisation)

    def __str__(self):
        if self.normalised == "user":
            return f"Broken power law IMF with transition mass {self.Mtransition}, and exponents {self.α1, self.α2}, normalised by {self.normalised} with constant {self.ξ0}."
        elif isinstance(self.normalised, tuple):
            if len(self.normalised) == 3 and self.normalised[0] in ("mass", "number"):
                return f"Broken power law IMF with with transition mass {self.Mtransition}, and exponents {self.α1, self.α2}, normalised by {self.normalised[0]} over the range {self.normalised[1:]} with constant {self.ξ0}."
            else:
                raise RuntimeError(f"Unexpected self.normalised flag in {self}")
        else:
            return f"with transition mass {self.Mtransition}, and exponents {self.α1, self.α2}. Not explicitly normalised."

    def __call__(self, M):
        """Calculate the value $\\xi(m)$."""
        # raise NotImplementedError
        if isinstance(M, (np.ndarray, tuple, list)):
            M = np.asarray(M)
            out = np.zeros_like(M, dtype=float)
            out[M <= self.Mtransition] = self.ξ0 * M[M <= self.Mtransition] ** self.α1
            out[M >= self.Mtransition] = self.ξ0 * self.ξcontinuity * M[M >= self.Mtransition] ** self.α2
            return out
        else:
            if M < self.Mtransition:
                return self.ξ0 * M ** self.α1
            else:
                return self.ξ0 * self.ξcontinuity * M ** self.α2

    def integrate(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
        $$
        \\int^{M_{max}}_{M_{min}} \\xi(m)dm =
        \\xi_0 \\frac{M_{max}^{\\alpha_1+1} - M_{min}^{\\alpha_1+1}}{\\alpha + 1} +
        \\xi_0 \\xi_\\textrm{continuity} \\frac{M_{max}^{\\alpha_2+1} - M_{min}^{\\alpha_2+1}}{\\alpha + 1},
        $$
        ```{note}
        The above equation doesn't take into account that Mmin and Mmax could both be above/below the transition mass, but the function call does.
        ```
        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        return self.ξ0 * (min(self.Mtransition, Mmax)**(self.α1+1) - min(self.Mtransition, Mmin)**(self.α1+1)) / (self.α1+1) \
             + self.ξ0 * self.ξcontinuity * (max(self.Mtransition, Mmax)**(self.α2+1) - max(self.Mtransition, Mmin)**(self.α2+1)) / (self.α2+1)

    def integrate_product(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
        $$
        \\int^{M_{max}}_{M_{min}} \\xi(m)dm =
        \\xi_0 \\frac{M_{max}^{\\alpha_1+2} - M_{min}^{\\alpha_1+2}}{\\alpha + 2} +
        \\xi_0 \\xi_\\textrm{continuity} \\frac{M_{max}^{\\alpha_2+2} - M_{min}^{\\alpha_2+2}}{\\alpha + 2},
        $$

        ```{note}
        The above equation doesn't take into account that Mmin and Mmax could both be above/below the transition mass, but the function call does.
        ```
        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        return self.ξ0 * (min(self.Mtransition, Mmax)**(self.α1+2) - min(self.Mtransition, Mmin)**(self.α1+2)) / (self.α1+2) \
             + self.ξ0 * self.ξcontinuity * (max(self.Mtransition, Mmax)**(self.α2+2) - max(self.Mtransition, Mmin)**(self.α2+2)) / (self.α2+2)

class L3IMF(InitialMassFunction):
    """
    Implement an L3 IMF of the form
    $$
    \\xi(m)dm = \\xi_0 \\left(\\frac{m}{\\mu}\\right)^{-\\alpha}\\left(1 + \\left(\\frac{m}{\\mu}\\right)^{1-\\alpha}\\right)^{-\\beta} dm,
    $$
    introduced in Maschberger (2013).

    You can specify if you want to normalise by mass or number, and the range.
    If a number is passed to normalisation then it is used
    If no normalisation is specified, then we take $\\xi_0 = 1$.

    Normalisation value provides the value we are normalising to, i.e. if you want the IMF to represent 100 solar masses.

    With the default parameter values of 2.3, 1.4, 0.2 this represents a best fit to Kroupa/Charier IMFs, provided in Maschberger (2013).
    """
    def __init__(self, alpha=2.3, beta=1.4, mu=0.2, normalisation=None, normalisation_value=1, Mmin=0.1, Mmax=100):
        self.α = alpha
        self.β = beta
        self.μ = mu

        # The constants end up a bit different to the Maschberger paper,
        # and I create a closure for B, which is similar to the incomplete beta function.
        a = (2 - self.α) / (1 - self.α)
        b = self.β - a
        self.integral_constant = self.μ / (1 - self.α) / (1 - self.β)
        self.integral_product_constant = self.μ**2 / (1 - self.α) * special.beta(a, b)
        self._B = self._create_B_closure(a, b)

        self.ξ0 = 1
        self.normalised = False
        self.Mmin = Mmin
        self.Mmax = Mmax

        if normalisation == "mass":
            self.normalise_by_mass(Mmin, Mmax, normalisation_value)
        elif normalisation == "number":
            self.normalise_by_number(Mmin, Mmax, normalisation_value)
        elif isinstance(normalisation, (int, float)):
            self.set_normalisation(normalisation)

    def __str__(self):
        if self.normalised == "user":
            return f"L3 IMF with parameters {self.α, self.β, self.μ}, normalised by {self.normalised} with constant {self.ξ0}."
        elif isinstance(self.normalised, tuple):
            if len(self.normalised) == 3 and self.normalised[0] in ("mass", "number"):
                return f"L3 IMF with parameters {self.α, self.β, self.μ}, normalised by {self.normalised[0]} over the range {self.normalised[1:]} with constant {self.ξ0}."
            else:
                raise RuntimeError(f"Unexpected self.normalised flag in {self}")
        else:
            return f"L3 IMF with parameters {self.α, self.β, self.μ}. Not explicitly normalised."

    def _create_B_closure(self, a, b):
        """
        Create a closure for the beta function as in Maschberger (2013), so we don't need to carry around a, b params.
        """
        def _B(M):
            """B integral (almost) as in Maschberger (2013)."""
            t = (M / self.μ)**(1-self.α) / (1 + (M / self.μ)**(1-self.α))
            return special.betainc(a, b, t)

        return _B

    def _G(self, M):
        """
        Auxilliary function from Maschberger (2013) used in other calculations.

        Has the form
        $$
        G(m) = \\left[1 + \\left(\\frac{m}{\\mu}\\right)^{1-\\alpha} \\right]^{1-\\beta}
        $$
        """
        return (1 + (M / self.μ)**(1-self.α))**(1 - self.β)

    def __call__(self, M):
        """Calculate the value $\\xi(m)$."""
        return self.ξ0 * (M / self.μ)**-self.α * (1 + (M / self.μ)**(1-self.α))**(-self.β)

    def integrate(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
            $\\int^{M_{max}}_{M_{min}} \\xi(m)dm = .

        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        return self.ξ0 * self.integral_constant * (self._G(Mmax) - self._G(Mmin))

    def integrate_product(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
            $\\int^{M_{max}}_{M_{min}} m\\xi(m)dm = ....

        From equation 25 of Maschberger (2013).

        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        return self.ξ0 * self.integral_product_constant * (self._B(Mmax) - self._B(Mmin))

#####################################################################################################
# Population III IMFs                                                                               #
# Power law IMFs are also commonly used for population III stars, but with different slopes/limits. #
# Some people use lognormal functions, similar to a Chabrier IMF without the power law tail.        #
# It is also possible to have a power law slope with some expontential cutoff.                      #
#####################################################################################################
# Like that in e.g. Larson (1973); Tumlinson (2006)
class LognormalIMF(InitialMassFunction):
    """
    Implement a Lognormal IMF of the form
    $$
    \\xi(m)dm = \\xi_0\\frac{1}{m}\\exp\\left(\\frac{-\\left(\\log m - \\log m_c\\right)^2}{2\\sigma^2}\\right)
    $$

    You can specify if you want to normalise by mass or number, and the range.
    If a number is passed to normalisation then it is used
    If no normalisation is specified, then we take $\\xi_0 = 1$.

    Normalisation value provides the value we are normalising to, i.e. if you want the IMF to represent 100 solar masses.

    Default values are not to be trusted blindly when considering Population III! The IMF is highly uncertain and you will need to make a decision.
    """
    def __init__(self, mc, sigma=1, normalisation=None, normalisation_value=1, Mmin=1, Mmax=500):
        self.mc = mc
        self.σ = sigma
        # Constants that go in front of integrals.
        self.lognormal_integral_constant = self.σ * np.log(10) * np.sqrt(np.pi/2)
        self.lognormal_product_integral_constant = np.sqrt(np.pi/2) * self.mc * self.σ * np.exp(np.log(10)**2*0.5*self.σ**2) * np.log(10)

        self.ξ0 = 1
        self.normalised = False
        self.Mmin = Mmin
        self.Mmax = Mmax

        if normalisation == "mass":
            self.normalise_by_mass(Mmin, Mmax, normalisation_value)
        elif normalisation == "number":
            self.normalise_by_number(Mmin, Mmax, normalisation_value)
        elif isinstance(normalisation, (int, float)):
            self.set_normalisation(normalisation)

    def __str__(self):
        if self.normalised == "user":
            return f"Lognormal IMF with characteristic mass {self.mc}, std. dev. {self.σ}, normalised by {self.normalised} with constant {self.ξ0}."
        elif isinstance(self.normalised, tuple):
            if len(self.normalised) == 3 and self.normalised[0] in ("mass", "number"):
                return f"Lognormal IMF with characteristic mass {self.mc}, std. dev. {self.σ}, normalised by {self.normalised[0]} over the range {self.normalised[1:]} with constant {self.ξ0}."
            else:
                raise RuntimeError(f"Unexpected self.normalised flag in {self}")
        else:
            return f"Lognormal IMF with characteristic mass {self.mc}, std. dev. {self.σ}. Not explicitly normalised."

    def __call__(self, M):
        """Calculate the value $\\xi(m)$."""
        # raise NotImplementedError
        # Vectorised as it is
        return self.ξ0 * np.exp( - (np.log10(M / self.mc))**2 / (2*self.σ**2)) / M

    def integrate(self, Mmin, Mmax):
        """
        Gives the analytical[1] solution to the integral:
            $\\int^{M_{max}}_{M_{min}} \\xi(m)dm = FILL IN$.

        NOTE: For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.

        Notes
        -----
        [1] There is no true analytical solution to the integral of the log-normal function, instead being represented in terms of the error function.
        """
        # NOTE: This mu is different to the mu in self.integrate_product
        μ1 = np.log10(Mmin / self.mc) / (np.sqrt(2)*self.σ)
        μ2 = np.log10(Mmax / self.mc) / (np.sqrt(2)*self.σ)
        return self.ξ0 * self.lognormal_integral_constant * (special.erf(μ2) - special.erf(μ1))

    def integrate_product(self, Mmin, Mmax):
        """
        Gives the analytical[1] solution to the integral:
            $\\int^{M_{max}}_{M_{min}} m\\xi(m)dm = FILL IN$.

        NOTE: For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.

        Notes
        -----
        [1] There is no true analytical solution to the integral of the log-normal function, instead being represented in terms of the error function.
        """
        # This took a while and lots of wrong factors and stuff, so please refer to computer algebra immediately: https://www.integral-calculator.com/#expr=exp%28-%28log_10%28x%29%20-%20log_10%28c%29%29%5E2%2F%282sigma%5E2%29&simplify=1
        # NOTE: This mu is different to the mu in self.integrate
        μ1 = (np.log10(Mmin/self.mc) - self.σ**2*np.log(10)) / (np.sqrt(2)*self.σ)
        μ2 = (np.log10(Mmax/self.mc) - self.σ**2*np.log(10)) / (np.sqrt(2)*self.σ)
        return self.ξ0 * self.lognormal_product_integral_constant * (special.erf(μ2) - special.erf(μ1))

# Like that in Wise+ (2012)
class GeneralisedGammaIMF(InitialMassFunction):
    """
    Implement a Generalised Gamma Function IMF of the form $\\xi(m) = \\xi_0 m^\\alpha\\exp\\left[-\\left(\\frac{m}{m_c}\\right)^\\beta\\right]$.

    You can specify if you want to normalise by mass or number, and the range.
    If a number is passed to normalisation then it is used
    If no normalisation is specified, then we take $\\xi_0 = 1$.

    Normalisation value provides the value we are normalising to, i.e. if you want the IMF to represent 100 solar masses.

    The default values of alpha, beta, mc come from [Wise+ (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...745...50W/abstract).
    """
    def __init__(self, alpha=-2.3, beta=-1.6, mc=10, normalisation=None, normalisation_value=1, Mmin=0.1, Mmax=100):
        self.α = alpha
        self.β = beta
        self.mc = mc

        self.ξ0 = 1
        self.normalised = False
        self.Mmin = Mmin
        self.Mmax = Mmax

        # Useful to define for the integrals - is it?
        # Needs to be defined before the normalisation step too
        self.integral_constant = self.mc**(self.α + 1) / self.β

        if normalisation == "mass":
            self.normalise_by_mass(Mmin, Mmax, normalisation_value)
        elif normalisation == "number":
            self.normalise_by_number(Mmin, Mmax, normalisation_value)
        elif isinstance(normalisation, (int, float)):
            self.set_normalisation(normalisation)

    @staticmethod
    def Γ_lower_incomplete(n_plus_1, u):
        """
        Lower incomplete gamma function for a given exponent.

        Defined as $\\Gamma(n+1, u) = \\int_0^u t^n e^{-t}dt$.
        """
        # Scipy provides a regularised lower incomplete gamma function, so the intergral we want, divided by the normal gamma function
        return special.gamma(n_plus_1) * special.gammainc(n_plus_1, u)

    def __str__(self):
        if self.normalised == "user":
            return f"Exponentially Truncated Power Law IMF with exponent {self.α}, characteristic mass {self.mc}, and damping exponent {self.β}, normalised by {self.normalised} with constant {self.ξ0}."
        elif isinstance(self.normalised, tuple):
            if len(self.normalised) == 3 and self.normalised[0] in ("mass", "number"):
                return f"Exponentially Truncated Power Law IMF with exponent {self.α}, characteristic mass {self.mc}, and damping exponent {self.β}, normalised by {self.normalised[0]} over the range {self.normalised[1:]} with constant {self.ξ0}."
            else:
                raise RuntimeError(f"Unexpected self.normalised flag in {self}")
        else:
            return f"Exponentially Truncated Power Law IMF with exponent {self.α}, characteristic mass {self.mc}, and damping exponent {self.β}. Not explicitly normalised."

    def __call__(self, M):
        """Calculate the value $\\xi(m) = \\xi_0 m^\\alpha\\exp\\left[-\\left(\\frac{m}{m_c}\\right)^\\beta\\right]$."""
        return self.ξ0 * M**self.α * np.exp(-(M / self.mc)**self.β)

    def integrate(self, Mmin, Mmax):
        """
        Gives the analytical solution to the integral:
            $\\int^{M_{max}}_{M_{min}} \\xi(m)dm = \\xi_0 \\frac{M_c^{\\alpha+1}}{\\beta}\\left[\\Gamma(n+1, u_{max}) - \\Gamma(n+1, u_{min})\\right]$.

        Where $\\Gamma(n+1, u) = \\int_0^u t^n e^{-t}dt$ is the lower incomplete gamma function, $n+1 = \\frac{\\alpha+1}{\\beta}$, and $u_X = \\left(\\frac{M_X}{M_c}\\right)^{\\beta}$.

        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        u1 = (Mmin / self.mc) ** self.β
        u2 = (Mmax / self.mc) ** self.β
        n_plus_1 = (self.α + 1) / self.β
        return self.ξ0 * self.integral_constant * (self.Γ_lower_incomplete(n_plus_1, u2) - self.Γ_lower_incomplete(n_plus_1, u1))

    def integrate_product(self, Mmin, Mmax):
        """
       Gives the analytical solution to the integral:
            $\\int^{M_{max}}_{M_{min}} \\xi(m)dm = \\xi_0 \\frac{M_c^{\\alpha+2}}{\\beta}\\left[\\Gamma(n+1, u_{max}) - \\Gamma(n+1, u_{min})\\right]$.

        Where $\\Gamma(n+1, u) = \\int_0^u t^n e^{-t}dt$ is the lower incomplete gamma function, $n+1 = \\frac{\\alpha+2}{\\beta}$, and $u_X = \\left(\\frac{M_X}{M_c}\\right)^{\\beta}$.

        ```{warning}
        For performance reasons we do not check Mmin < Mmax, or are within the bounds provided for normalisation.
        ```
        """
        u1 = (Mmin / self.mc) ** self.β
        u2 = (Mmax / self.mc) ** self.β
        n_plus_1 = (self.α + 2) / self.β
        return self.ξ0 * self.mc * self.integral_constant * (self.Γ_lower_incomplete(n_plus_1, u2) - self.Γ_lower_incomplete(n_plus_1, u1))

class ExponentialCutoffPowerLawIMF(GeneralisedGammaIMF):
    def __init__(self, *args, **kwargs):
        raise DeprecationWarning("`ExponentialCutoffPowerLawIMF` has been renamed `GeneralisedGammaIMF`.")

def popII_IMFs():
    import matplotlib.pyplot as plt
    from scipy.integrate import quad
    M = np.geomspace(0.1, 100, 100)
    for imf in [
        PowerLawIMF(normalisation="mass"),
        BrokenPowerLawIMF(normalisation="mass"),
        ChabrierIMF(normalisation="mass"),
        L3IMF(normalisation="mass")
            ]:
        print(imf)
        print("Difference between true mass and integrated (closer to 0 better): ", imf.integrate_product(imf.Mmin, imf.Mmax) - 1)

        fig, axes = plt.subplot_mosaic([["IMF value", "IMF value"],["Integrated IMF","Integrated IMF - scipy quad"],["Integrated m*IMF","Integrated m*IMF - scipy quad"]], sharex=True)
        axes["IMF value"].loglog(M, imf(M))
        axes["IMF value"].set(xlabel="Mass", ylabel="IMF value", title="IMF")

        axes["Integrated IMF"].loglog(M[1:], [imf.integrate(M[0], MM) for MM in M[1:]])
        axes["Integrated IMF"].set(xlabel="Mass", ylabel="Integrated IMF", title="Integrated IMF")

        axes["Integrated IMF - scipy quad"].loglog(M[1:], [abs(imf.integrate(M[0], MM) - quad(imf, M[0], MM, limit=100000000, epsabs=1e-14, epsrel=1e-14)[0]) for MM in M[1:]])
        axes["Integrated IMF - scipy quad"].set(xlabel="Mass", ylabel="abs(difference)", title="Integrated IMF - scipy quad")

        axes["Integrated m*IMF"].loglog(M[1:], [imf.integrate_product(M[0], MM) for MM in M[1:]])
        axes["Integrated m*IMF"].set(xlabel="Mass", ylabel="Integrated m*IMF", title="Integrated m*IMF")

        axes["Integrated m*IMF - scipy quad"].loglog(M[1:], [abs(imf.integrate_product(M[0], MM) - quad(lambda m: m*imf(m), M[0], MM, limit=100000000, epsabs=1e-14, epsrel=1e-14)[0]) for MM in M[1:]])
        axes["Integrated m*IMF - scipy quad"].set(xlabel="Mass", ylabel="abs(difference)", title="Integrated m*IMF - scipy quad")

        plt.suptitle(imf)
        plt.show()

    for imf in [PowerLawIMF(normalisation="mass"),
                BrokenPowerLawIMF(normalisation="mass"),
                ChabrierIMF(normalisation="mass"),
                L3IMF(normalisation="mass")]:
        plt.loglog(M, imf(M))
    plt.title("Number count of stars")
    plt.show()
    for imf in [PowerLawIMF(normalisation="mass"),
                BrokenPowerLawIMF(normalisation="mass"),
                ChabrierIMF(normalisation="mass"),
                L3IMF(normalisation="mass")]:
        plt.loglog(M[1:], [imf.integrate_product(M[0], MM) for MM in M[1:]])
    plt.title("Cumulative Mass")
    plt.show()

def popIII_IMFs():
    import matplotlib.pyplot as plt
    from scipy.integrate import quad
    M = np.geomspace(1, 500, 100)
    for imf in [
        # PowerLawIMF(normalisation="mass", Mmin=1, Mmax=100),
        PowerLawIMF(normalisation="mass", Mmin=10, Mmax=500),
        # LognormalIMF(normalisation="mass", Mmin=1, Mmax=100, mc=10),
        LognormalIMF(normalisation="mass", Mmin=1, Mmax=500, mc=30),
        GeneralisedGammaIMF(normalisation="mass", Mmin=1, Mmax=100, mc=10)
            ]:
        print(imf)
        print("Difference between true mass and integrated (closer to 0 better): ", imf.integrate_product(imf.Mmin, imf.Mmax) - 1)

        fig, axes = plt.subplot_mosaic([["IMF value", "IMF value"],["Integrated IMF","Integrated IMF - scipy quad"],["Integrated m*IMF","Integrated m*IMF - scipy quad"]], sharex=True)
        axes["IMF value"].loglog(M, imf(M))
        axes["IMF value"].set(xlabel="Mass", ylabel="IMF value", title="IMF")

        axes["Integrated IMF"].loglog(M[1:], [imf.integrate(M[0], MM) for MM in M[1:]])
        axes["Integrated IMF"].set(xlabel="Mass", ylabel="Integrated IMF", title="Integrated IMF")

        axes["Integrated IMF - scipy quad"].loglog(M[1:], [abs(imf.integrate(M[0], MM) - quad(imf, M[0], MM, limit=100000000, epsabs=1e-14, epsrel=1e-14)[0]) for MM in M[1:]])
        axes["Integrated IMF - scipy quad"].set(xlabel="Mass", ylabel="abs(difference)", title="Integrated IMF - scipy quad")

        axes["Integrated m*IMF"].loglog(M[1:], [imf.integrate_product(M[0], MM) for MM in M[1:]])
        axes["Integrated m*IMF"].set(xlabel="Mass", ylabel="Integrated m*IMF", title="Integrated m*IMF")

        axes["Integrated m*IMF - scipy quad"].loglog(M[1:], [abs(imf.integrate_product(M[0], MM) - quad(lambda m: m*imf(m), M[0], MM, limit=100000000, epsabs=1e-14, epsrel=1e-14)[0]) for MM in M[1:]])
        axes["Integrated m*IMF - scipy quad"].set(xlabel="Mass", ylabel="abs(difference)", title="Integrated m*IMF - scipy quad")

        plt.suptitle(imf)
        plt.show()

    plot_IMFs = [
        PowerLawIMF(normalisation="mass", Mmin=1, Mmax=100),
        PowerLawIMF(normalisation="mass", Mmin=10, Mmax=500),
        LognormalIMF(normalisation="mass", Mmin=1, Mmax=100, mc=10),
        LognormalIMF(normalisation="mass", Mmin=1, Mmax=500, mc=30),
        GeneralisedGammaIMF(normalisation="mass", Mmin=1, Mmax=100, mc=10)
            ]

    for imf in plot_IMFs:
        mask = (imf.Mmin <= M) & (M <= imf.Mmax)
        plt.loglog(M[mask], imf(M[mask]), label=imf)
    plt.legend()
    plt.title("Number count of stars")
    plt.show()
    for imf in plot_IMFs:
        plt.loglog(M[1:], [imf.integrate_product(max(M[0], imf.Mmin), MM) if MM <= imf.Mmax else None for MM in M[1:]])
    plt.title("Cumulative Mass")
    plt.show()

def main():
    popII_IMFs()
    popIII_IMFs()

if __name__ == "__main__":
    main()