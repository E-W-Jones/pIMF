# The Theory

This part of the documentation is about the theory behind the IMF and choices about implementation.

```{contents}

```

## Basics

The stellar initial mass function (IMF) describes the distribution of stellar masses that form. It has wide-ranging effects and is currently one of the largest uncertainties in our galaxy models, due to how hard it is to determine it reliably. From observations one has to detangle the effects of stellar lifetimes and the star formation history of a galaxy. In simulations there is an inherent trade off between the number of stars you can form, the smallest masses resolved, and the timescale that can be studied.

While any definition of the IMF represents the same basic idea, the details introduce subtle changes, adding or removing factors of mass. In this package we use the number of stars per (linear) mass bin $$\xi(M) = \frac{dN}{dM}.$$

It is also popular to define the IMF based on the mass of stars per mass bin: $$\Xi(M) = M\xi(M),$$ or using the logarithm of mass: $$\Xi(M) = \frac{dN}{d\log M} = M\xi(M),$$ where we have demonstrated the equivalence of different definitions using the chain rule and the fact $\frac{dM}{d\log M} = M$.

### Why is the IMF important?

The IMF is important as it determines how many stars form at a given mass. This is a necessary assumption to make for any method of linking observations, where individual stars cannot be resolved,  to models, which start with single stars at a given mass. And, since the initial mass of a star has a large effect on its evolution---including how bright it is, how many ionising photons it emits, how long it lives, and whether or not it ends its life as a supernova---a small change in the number of massive stars can cause a large change on the properties of a galaxies.

### Further Reading

There is a wealth of resources that I can't even begin to cover. Two recent reviews are [Kroupa+ (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv241007311K/abstract), a textbook chapter that covers some of the theory and observational evidence; and [Hennebelle & Grudić (2024)](https://ui.adsabs.harvard.edu/abs/2024ARA%26A..62...63H/abstract), which covers more of the simulation side.

### Population II vs Population III
Population II stars are stars that form with metals. Population III stars are stars that form without metals. ([This notation is historical](https://en.wikipedia.org/wiki/Stellar_population#Chemical_classification_by_Walter_Baade), and therefore not the most intuitive.)

The population II IMF is thought to be relatively well constrained. The two most commonly used parameterisations are that of [Kroupa (2001)](functional_forms/kroupa.md) and [Chabrier (2003)](functional_forms/chabrier.md). These are preferred to the [Salpeter (1955)](functional_forms/salpeter.md) model as they turnover at low masses, while having the same power law at the high-mass end. However, the Salpeter IMF is still used when studying the most massive stars.

Population III stars form out of metal-free gas, compared to population II which are enriched in some way. How this, alongside other conditions of the early universe, changes the IMF is a currently open area of astronomy research. It is complicated by the fact that, alongside everything mentioned before, population III stars have never been observed (and possibly never will). This means everything we know about them must be gleaned from their impact on their surroundings, like the abundance patterns of elements they produce.

In this package we include a variety of IMFs that are suitable for population II or population III stars, and are always happy to include more.

## Integrals and Observables
The most common way to use the IMF is to treat it as a probablity density function, and use it to weight averages of functions of mass.

### Simple Integral
The number of stars between $M_\mathrm{min}$, $M_\mathrm{max}$ is given by $$\int_{M_\mathrm{min}}^{M_\mathrm{max}}\xi(M)\mathrm{d}M.$$

```{admonition} Example: Number of Supernovae
:class: toggle, tip
To calculate the number of supernovae that go off, assuming a [Salpeter (1955)](./functional_forms/salpeter.md) IMF, integrate from $8M_\odot$ to $100M_\odot$:
\begin{align*}
N_\mathrm{SN} &= \int_{8M_\odot}^{100M_\odot}\xi_0 M^{-2.35}\mathrm{d}M\\
&= \xi_0\frac{100M_\odot^{-1.35} - 8M_\odot^{-1.35}}{-1.35}\\
&\approx0.04M_\odot^{-1.35}\xi_0
\end{align*}
```

This integral can be calculated with the `integrate` method of your IMF of choice.

### Mass-weighted Integral
The mass in stars between $M_\mathrm{min}$, $M_\mathrm{max}$ is given by $$\int_{M_\mathrm{min}}^{M_\mathrm{max}}M\xi(M)\mathrm{d}M.$$

This integral can be calculated with the `integrate_product` method of your IMF of choice.

### Arbitrary Functions
Consider an arbitrary quantity you are interested in that is a function of stellar mass. For example, the luminosity $L(M)$. 

This means the total luminosity will be given by the luminosity of a star at a given mass, multiplied by the number of stars that have that mass: $$\int^{M_\textrm{max}}_{M_\textrm{min}} L(M)\xi(M)dM.$$

Alternatively, you might be interested in the average luminosity, which is likely as you don't have a specific number (or mass) of stars in mind. This would be given by
$$\frac{\int^{M_\textrm{max}}_{M_\textrm{min}} L(M)\xi(M)dM}{\int^{M_\textrm{max}}_{M_\textrm{min}} \xi(M)dM}.$$

## Functional Forms
We currently have several functional forms interested. If there's one you'd like to see included, check out our [subclassing guide](subclassing.ipynb), and consider submitting a pull request to get it included! The below pages include the functional form, some references, notes on implementation, and a derivation of the [simple](#simple-integral) and [mass weighted integrals](#mass-weighted-integral) implemented in the code.
```{toctree}
:maxdepth: 2
functional_forms/index.md
```