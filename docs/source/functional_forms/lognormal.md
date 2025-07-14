# Lognormal
## Form
Suggested by, e.g., [Larson (1973)](https://ui.adsabs.harvard.edu/abs/1973MNRAS.161..133L/abstract), the IMF takes the form of a [lognormal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution):
$$\xi(M) = \xi_0\frac{1}{M}\exp\left(\frac{-\left(\log M - \log M_c\right)^2}{2\sigma^2}\right).$$

## Implemented in
This is implemented in the {py:class}`~pimf.initialmassfunction.LognormalIMF` class, where the user is free to pick whatever parameters they like or use the default (same as above).

## Integrals
Behind the scenes we will use [`scipy.special.erf`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf.html) to evaluate these integrals. Defining $$\mathtt{erf}(z) = \frac{2}{\sqrt{\pi}} \int^{z}_{0} e^{-t^2}\mathrm{d}t$$ and recalling $\int^{b}_{a}f(x)\mathrm{d}x = \int^{b}_{0}f(x)\mathrm{d}x - \int^{a}_{0}f(x)\mathrm{d}x$ we can write
$$\begin{align*}
\int^{b}_{a} e^{-t^2}\mathrm{t}
    &= \int^{b}_{0} e^{-t^2}\mathrm{t} - \int^{a}_{0} e^{-t^2}\mathrm{t} \\
    &= \frac{\sqrt{\pi}}{2}\left(\mathtt{erf}(b) - \mathtt{erf}(a)\right) \\
\end{align*}$$


### Total number of stars
Use the substitution $$\mu = \frac{\log M - \log M_c}{\sqrt{2}\sigma} = \frac{\ln(M/M_c)}{\sqrt{2}\sigma\ln10}$$ which has the derivative $$\mathrm{d}\mu = \frac{1}{\sqrt{2}\sigma\ln10}\frac{\mathrm{d}M}{M}$$ to change variables $$\xi(M)\mathrm{d}M = \xi_0\sqrt{2}\sigma\ln10\exp\left(-\mu^2\right)\mathrm{d}\mu$$

$$\begin{align*}
\int^{M_\textrm{max}}_{M_\textrm{min}} \xi(M)\mathrm{d}M
    &= \int^{\mu_\textrm{max}}_{\mu_\textrm{min}} \xi(\mu)\mathrm{d}\mu \\
    &= \xi_0\sqrt{2}\sigma\ln10\int^{\mu_\textrm{max}}_{\mu_\textrm{min}} \exp\left(-\mu^2\right)\mathrm{d}\mu \\
    &= \xi_0\sigma\ln10\sqrt{\frac{\pi}{2}}\left(\mathtt{erf}(\mu_\textrm{max}) - \mathtt{erf}(\mu_\textrm{min})\right)
\end{align*}$$

### Total mass of stars
We will have an extra factor of $M$. Rearranging our substitution to $M = M_c e^{\sqrt{2}\sigma\ln 10\mu}$:
$$\begin{align*}
\int^{M_\textrm{max}}_{M_\textrm{min}} M\xi(M)\mathrm{d}M
    &= \int^{\mu_\textrm{max}}_{\mu_\textrm{min}} M_c e^{\sqrt{2}\sigma\ln 10\mu}\xi(\mu)\mathrm{d}\mu \\
    &= \xi_0\sqrt{2}\sigma\ln10M_c\int^{\mu_\textrm{max}}_{\mu_\textrm{min}} e^{\sqrt{2}\sigma\ln 10\mu}\exp\left(-\mu^2\right)\mathrm{d}\mu \\
    &= \xi_0\sqrt{2}\sigma\ln10M_c\int^{\mu_\textrm{max}}_{\mu_\textrm{min}} \exp\left(-\left[\mu^2 - \sqrt{2}\sigma\ln 10\mu\right]\right)\mathrm{d}\mu 
\end{align*}$$

This can be solved by another substitution, which we use to complete the square
$$\begin{align*}
\mu' &= \mu - \frac{\sigma\ln10}{\sqrt{2}} \\
\Rightarrow \mu'^2 -\frac{(\sigma\ln10)^2}{2} &= \mu^2 - \sqrt{2}\sigma\ln10\mu \\
\int^{M_\textrm{max}}_{M_\textrm{min}} M\xi(M)\mathrm{d}M
    &= \xi_0\sqrt{2}\sigma\ln10M_c\int^{\mu_\textrm{max}}_{\mu_\textrm{min}} \exp\left(-\left[\mu'^2 -\frac{(\sigma\ln10)^2}{2}\right]\right)\mathrm{d}\mu \\
    &= \xi_0\sqrt{2}\sigma\ln10M_c e^\frac{(\sigma\ln10)^2}{2}\int^{\mu_\textrm{max}}_{\mu_\textrm{min}} \exp\left(-\mu'^2\right)\mathrm{d}\mu \\
    &= \xi_0\sigma\ln10M_c e^\frac{(\sigma\ln10)^2}{2}\sqrt{\frac{\pi}{2}}\left(\mathtt{erf}(\mu'_\textrm{max}) - \mathtt{erf}(\mu'_\textrm{min})\right)
\end{align*}$$


:::{seealso}
* [Chabrier (2003)](./chabrier.md)
:::