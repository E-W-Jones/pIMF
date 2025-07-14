# Power Law (Salpeter 1955)
## Form
[Salpeter (1955)](https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract) suggest a simple power law: $$\xi(M) = \xi_0M^\alpha$$ with $$\alpha = -2.35.$$

## Implemented in
This is implemented in the {py:class}`~pimf.initialmassfunction.PowerLawIMF` class, where the user is free to pick whatever slope they like or use the default (same as above).

## Integrals

### Total number of stars
$$\begin{align}
\int^{M_\textrm{max}}_{M_\textrm{min}}\xi(M)\mathrm{d}M 
    &= \xi_0\int^{M_\textrm{max}}_{M_\textrm{min}}M^\alpha\mathrm{d}M \\
    &= \xi_0\left[\frac{M^{\alpha+1}}{\alpha + 1}\right]^{M_\textrm{max}}_{M_\textrm{min}} \\
    &= \xi_0\frac{M^{\alpha+1}_\textrm{max} - M^{\alpha+1}_\textrm{min}}{\alpha + 1}.
\end{align}$$

### Total mass of stars
$$\begin{align}
\int^{M_\textrm{max}}_{M_\textrm{min}}M\xi(M)\mathrm{d}M 
    &= \xi_0\int^{M_\textrm{max}}_{M_\textrm{min}}M^{\alpha+1}\mathrm{d}M \\
    &= \xi_0\left[\frac{M^{\alpha+2}}{\alpha + 2}\right]^{M_\textrm{max}}_{M_\textrm{min}} \\
    &= \xi_0\frac{M^{\alpha+2}_\textrm{max} - M^{\alpha+2}_\textrm{min}}{\alpha + 2}.
\end{align}$$

These are only valid for $\alpha \neq -1$ and $\alpha \neq -2$, respectively. In either of these cases the integrals evaluate to $$\xi_0\ln\left(\frac{M_\textrm{max}}{M_\textrm{min}}\right).$$

:::{seealso}
* [Broken Power Law (Kroupa 2001)](./kroupa.md)
* [Chabrier (2003)](./chabrier.md)
:::