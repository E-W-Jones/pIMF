# Generalised Gamma (Wise+ 2012)
```{danger}
This IMF used to be called `ExponentialCutoffPowerLawIMF`. It's possible some old references to this still exist.
```
## Form
Also referred to as the _generalized
Rosin-Rammler function_ by [Chabrier (2003)](https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract), the [Generalised Gamma Distribution](https://en.wikipedia.org/wiki/Generalized_gamma_distribution) can be expressed as $$\xi(m) = \xi_0 M^\alpha\exp\left[-\left(\frac{M}{M_c}\right)^\beta\right].$$ This looks like a power law with an exponential cutoff at the low-mass end, which makes it tempting for studies of [Population III stars](../theory.md#population-ii-vs-population-iii).

The default values in our implementation are $$\alpha = -2.3, \beta=1.6, M_c = 10M_\odot.$$

[Wise+ (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...745...50W/abstract) provide the values $$\alpha = -2.3, \beta=1.6,$$ but suggested a characteristic mass $M_c\sim100M_\odot$. However, recent studies take values closer to $10M_\odot$ (e.g. [Brauer+ (2025)](https://ui.adsabs.harvard.edu/abs/2025ApJ...980...41B/abstract)). 

## Implemented in
This is implemented in the {py:class}`~pimf.initialmassfunction.GeneralisedGammaIMF` class, where the user is free to pick whatever parameters they like or use the default (same as above).

## Integrals
Under the hood we will be using `scipy`'s [`gamma`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gamma.html) and [`gammainc`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammainc.html) to compute the lower incomplete gamma function $$\mathtt{gamma}(n+1)\mathtt{gammainc}(n+1, u) = \Gamma(n+1, u) = \int_0^u t^n e^{-t}\mathrm{d}t.$$

We can change the lower limits using the fact $\int^{b}_{a}f(x)\mathrm{d}x = \int^{b}_{0}f(x)\mathrm{d}x - \int^{a}_{0}f(x)\mathrm{d}x$:
$$
\int_\ell^u t^n e^{-t}\mathrm{d}t = \Gamma(n+1, u) - \Gamma(n+1, \ell)
$$

We will also use the substitution $$u = \left(\frac{M}{M_c}\right)^\beta \Rightarrow\mathrm{d}M = \frac{M_c}{\beta}u^{\frac{1-\beta}{\beta}}\mathrm{d}u.$$
### Total number of stars
$$\begin{align}
\int^{M_\textrm{max}}_{M_\textrm{min}}\xi(M)\mathrm{d}M 
    &= \xi_0\int^{M_\textrm{max}}_{M_\textrm{min}}M^\alpha\exp\left\{-\left(\frac{M}{M_c}\right)^\beta\right\}\mathrm{d}M \\
    &= \frac{\xi_0M_c^{\alpha+1}}{\beta}\int^{u(M_\textrm{max})}_{u(M_\textrm{min})}u^{\frac{\alpha+1}{\beta}-1}e^{-u}\mathrm{d}u \\
    &= \frac{\xi_0M_c^{\alpha+1}}{\beta}\left\{ \Gamma\left(\frac{\alpha+1}{\beta}, \left[\frac{M_\mathrm{max}}{M_c}\right]^\beta\right) - \Gamma\left(\frac{\alpha+1}{\beta}, \left[\frac{M_\mathrm{min}}{M_c}\right]^\beta\right) \right\}.
\end{align}$$

### Total mass of stars
$$\begin{align}
\int^{M_\textrm{max}}_{M_\textrm{min}}\xi(M)\mathrm{d}M 
    &= \xi_0\int^{M_\textrm{max}}_{M_\textrm{min}}M^{\alpha+1}\exp\left\{-\left(\frac{M}{M_c}\right)^\beta\right\}\mathrm{d}M \\
    &= \frac{\xi_0M_c^{\alpha+2}}{\beta}\int^{u(M_\textrm{max})}_{u(M_\textrm{min})}u^{\frac{\alpha+2}{\beta}-1}e^{-u}\mathrm{d}u \\
    &= \frac{\xi_0M_c^{\alpha+2}}{\beta}\left\{ \Gamma\left(\frac{\alpha+2}{\beta}, \left[\frac{M_\mathrm{max}}{M_c}\right]^\beta\right) - \Gamma\left(\frac{\alpha+2}{\beta}, \left[\frac{M_\mathrm{min}}{M_c}\right]^\beta\right) \right\}.
\end{align}$$

:::{seealso}
* [Power Law (Salpeter 1955)](./salpeter.md)
:::