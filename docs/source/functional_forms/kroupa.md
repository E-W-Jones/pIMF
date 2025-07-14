# Broken Power Law (Kroupa 2001)
## Form
[Kroupa (2001)](https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract) suggested an extending the [Salpeter (1955)](./salpeter.md) IMF to a piecewise defined power law with several slopes:
$$
\xi(M) =
\xi_0\begin{cases}
M^{\alpha_1}, & M \leq M_\textrm{transition} \\
\xi_\textrm{continuity}M^{\alpha_2}, & M_\textrm{transition} \leq M
\end{cases}
$$
where for continuity one requires $$\xi_\textrm{continuity} = (M_\textrm{transition})^{\alpha_1 - \alpha_2}.$$

The values when fit to data are:
$$\begin{gather*}
\alpha_1 = -1.3 \\
\alpha_2 = -2.3 \\
M_\textrm{transisiton} = 0.5M_\odot
\end{gather*}$$


:::{note}
The original paper by Kroupa (2001) includes an additional mass range, below $0.08M_\odot$. This is below the limit where fusion can happen, so these objects are _brown dwarves_ instead of stars, and are ignored.
:::

## Implemented in
This is implemented in the {py:class}`~pimf.initialmassfunction.BrokenPowerLawIMF` class, where the user is free to pick whatever parameters they like or use the default (same as above).

## Integrals
As mentioned in [power law integrals](./salpeter.md#integrals), there are special case values of $\alpha$, that are not explicity mentioned here. I will also only consider the case $M_\textrm{min} < M_\textrm{transition} < M_\textrm{max}$. If this is not the case, simply ignore the term from that branch.

### Total number of stars
$$\begin{align*}
\int^{M_\textrm{max}}_{M_\textrm{min}} \xi(M)\mathrm{d}M
    &= \int^{M_\textrm{transition}}_{M_\textrm{min}} \xi_0 M^{\alpha_1}\textrm{d}M + \int^{M_\textrm{max}}_{M_\textrm{transition}} \xi_0\xi_\textrm{continuity}M^{\alpha_2} \textrm{d}M \\
    &= \left[\xi_0 \frac{M^{\alpha_1+1}}{\alpha_1+1}\right]^{M_\textrm{transition}}_{M_\textrm{min}} + \left[\xi_0\xi_\textrm{continuity}\frac{M^{\alpha_2+1}}{\alpha_2+1}\right]^{M_\textrm{max}}_{M_\textrm{transition}} \\
    &= \xi_0 \frac{M_\textrm{transition}^{\alpha_1+1} - M_\textrm{min}^{\alpha_1+1}}{\alpha_1 + 1} +
\xi_0 \xi_\textrm{continuity} \frac{M_\textrm{max}^{\alpha_2+1} - M_\textrm{transition}^{\alpha_2+1}}{\alpha_2 + 1}.
\end{align*}$$

### Total mass of stars
$$\begin{align*}
\int^{M_\textrm{max}}_{M_\textrm{min}} M\xi(M)\mathrm{d}M
    &= \int^{M_\textrm{transition}}_{M_\textrm{min}} \xi_0 M^{\alpha_1+1}\textrm{d}M + \int^{M_\textrm{max}}_{M_\textrm{transition}} \xi_0\xi_\textrm{continuity}M^{\alpha_2+1} \textrm{d}M \\
    &= \left[\xi_0 \frac{M^{\alpha_1+2}}{\alpha_1+2}\right]^{M_\textrm{transition}}_{M_\textrm{min}} + \left[\xi_0\xi_\textrm{continuity}\frac{M^{\alpha_2+2}}{\alpha_2+2}\right]^{M_\textrm{max}}_{M_\textrm{transition}} \\
    &= \xi_0 \frac{M_\textrm{transition}^{\alpha_1+2} - M_\textrm{min}^{\alpha_1+2}}{\alpha_1 + 2} +
\xi_0 \xi_\textrm{continuity} \frac{M_\textrm{max}^{\alpha_2+2} - M_\textrm{transition}^{\alpha_2+2}}{\alpha_2 + 2}.
\end{align*}$$

:::{seealso}
* [Power Law (Salpeter 1955)](./salpeter.md)
* [Chabrier (2003)](./chabrier.md)
:::
