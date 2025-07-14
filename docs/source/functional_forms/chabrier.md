# Chabrier (2003)
## Form
Similarly to [Kroupa (2001)](./kroupa.md), [Chabrier (2003)](https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract) extend the [Salpeter (1955)](./salpeter.md) power law. However, they introduce a lognormal distribution for the low-mass end, instead of power laws. This takes the form:

$$
\xi(M) = \xi_0
\begin{cases}
\frac{1}{M}\exp\left(\frac{-\left(\log M - \log M_c\right)^2}{2\sigma^2}\right), & M \leq 1M_\odot\\
\xi_\textrm{continuity} M^{\alpha}, & 1M_\odot \leq M
\end{cases}
$$
with best fit values
$$\begin{gather*}
M_c = 0.079\\
\sigma = 0.69 \\
\alpha = 2.3.
\end{gather*}$$

:::{note}
In theory one is free to choose any transition mass. However, this has not been implemented.
:::

For continuity, $$\xi_\textrm{continuity} = M^{-\alpha-1}_\textrm{transition}\exp\left(\frac{-\left(\log M_\textrm{transition} - \log M_c\right)^2}{2\sigma^2}\right).$$

## Implemented in
This is implemented in the {py:class}`~pimf.initialmassfunction.ChabrierIMF` class, where the user is free to pick whatever parameters they like or use the default (same as above).

## Integrals
As mentioned in [power law integrals](./salpeter.md#integrals), there are special case values of $\alpha$, that are not explicity mentioned here. I will also only consider the case $M_\textrm{min} < 1M_\odot < M_\textrm{max}$. If this is not the case, simply ignore the term from that branch.

For details of the lognormal integrals, including the substitutions, see [Lognormal integrals](./lognormal.md#integrals). We will use the `scipy` implementation of the error function, [`scipy.special.erf`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf.html): $$\mathtt{erf}(z) = \frac{2}{\sqrt{\pi}} \int^{z}_{0} e^{-t^2}\mathrm{d}t.$$

### Total number of stars
$$\begin{align*}
\int^{M_\textrm{max}}_{M_\textrm{min}} \xi(M)\mathrm{d}M
    &= \int^{1M_\odot}_{M_\textrm{min}} \frac{\xi_0}{M}\exp\left(\frac{-\left(\log M - \log M_c\right)^2}{2\sigma^2}\right)\mathrm{d}M
    + \int^{M_\textrm{max}}_{1M_\odot} \xi_0\xi_\textrm{continuity}M^{\alpha} \textrm{d}M \\
    &= \left[\xi_0 \sigma\ln10\sqrt{\frac{\pi}{2}}\mathtt{erf}(\mu)\right]^{\frac{\log(1M_\odot/M_c)}{\sqrt{2}\sigma}}_{\frac{\log(M_\textrm{min}/M_c)}{\sqrt{2}\sigma}}
    + \left[\xi_0\xi_\textrm{continuity}\frac{M^{\alpha+1}}{\alpha+1}\right]^{M_\textrm{max}}_{1M_\odot} \\
    &= \xi_0\sigma\ln10\sqrt{\frac{\pi}{2}}\left\{\mathtt{erf}\left(\frac{\log(1M_\odot/M_c)}{\sqrt{2}\sigma}\right) - \mathtt{erf}\left(\frac{\log(M_\textrm{min}/M_c)}{\sqrt{2}\sigma}\right)\right\} \\
    &\quad + \xi_0 \xi_\textrm{continuity} \frac{M_\textrm{max}^{\alpha+1} - 1M_\odot^{\alpha+1}}{\alpha + 1}.
\end{align*}$$

### Total mass of stars

$$\begin{align*}
\int^{M_\textrm{max}}_{M_\textrm{min}} M\xi(M)\mathrm{d}M 
    &= \int^{1M_\odot}_{M_\textrm{min}} \xi_0\exp\left(\frac{-\left(\log M - \log M_c\right)^2}{2\sigma^2}\right)\mathrm{d}M
    + \int^{M_\textrm{max}}_{1M_\odot} \xi_0\xi_\textrm{continuity}M^{\alpha+1} \textrm{d}M \\
    &= \left[\xi_0\sigma\ln10M_c e^\frac{(\sigma\ln10)^2}{2}\sqrt{\frac{\pi}{2}}\mathtt{erf}\left(\mu - \frac{\sigma\ln10}{\sqrt{2}}\right)\right]^{\frac{\log(1M_\odot/M_c)}{\sqrt{2}\sigma}}_{\frac{\log(M_\textrm{min}/M_c)}{\sqrt{2}\sigma}} \\
    &\quad + \left[\xi_0\xi_\textrm{continuity}\frac{M^{\alpha+2}}{\alpha+2}\right]^{M_\textrm{max}}_{1M_\odot} \\
    &= \xi_0\sigma\ln10\sqrt{\frac{\pi}{2}}\left\{\mathtt{erf}\left(\frac{\log(1M_\odot/M_c)}{\sqrt{2}\sigma} - \frac{\sigma\ln10}{\sqrt{2}}\right) - \mathtt{erf}\left(\frac{\log(M_\textrm{min}/M_c)}{\sqrt{2}\sigma} - \frac{\sigma\ln10}{\sqrt{2}}\right)\right\} \\
    &\quad + \xi_0 \xi_\textrm{continuity} \frac{M_\textrm{max}^{\alpha+2} - 1M_\odot^{\alpha+2}}{\alpha + 2}.
\end{align*}$$

:::{seealso}
* [Power Law (Salpeter 1955)](./salpeter.md)
* [Broken Power Law (Kroupa 2001)](./kroupa.md)
* [Lognormal](./lognormal.md)
:::
