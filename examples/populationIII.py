from pimf import PowerLawIMF, LognormalIMF, GeneralisedGammaIMF
import numpy as np
import matplotlib.pyplot as plt
from itertools import pairwise
# Package IMFs in namedtuples with info we want for plotting
from collections import namedtuple
imf = namedtuple("imf_instance", ["imf", "name", "plot_kwargs"])

def plot_helper(imf, num=50):
    M = np.geomspace(imf.Mmin, imf.Mmax, num=num)

TOTAL_MASS = 1

imfs = [
    imf(name="FiBY", imf=PowerLawIMF(normalisation="mass", normalisation_value=TOTAL_MASS, Mmin=21, Mmax=500), plot_kwargs=dict(color="C0")),
    imf(name="Sarmento+19", imf=LognormalIMF(60, normalisation="mass", normalisation_value=TOTAL_MASS), plot_kwargs=dict(color="C1")),
    imf(name="A-SLOTH", imf=PowerLawIMF(-1.77, normalisation="mass", normalisation_value=TOTAL_MASS, Mmin=13.6, Mmax=197), plot_kwargs=dict(color="C2")),
    imf(name="AEOS10", imf=GeneralisedGammaIMF(normalisation="mass", normalisation_value=TOTAL_MASS, Mmin=1, Mmax=100), plot_kwargs=dict(color="C3")),
    imf(name="AEOS20", imf=GeneralisedGammaIMF(normalisation="mass", normalisation_value=TOTAL_MASS, Mmin=1, Mmax=300, mc=20), plot_kwargs=dict(color="C3", linestyle="--")),
    imf(name="HIMALAYAS", imf=PowerLawIMF(normalisation="mass", normalisation_value=TOTAL_MASS, Mmin=10, Mmax=300), plot_kwargs=dict(color="C4"))
]

fig, axes = plt.subplots(2, 2)

for imf in imfs:
    Mmin, Mmax = imf.imf.Mmin, imf.imf.Mmax
    M = np.geomspace(Mmin, Mmax)
    # Number of stars
    axes[0, 0].plot(M, imf.imf(M), label=imf.name, **imf.plot_kwargs)
    axes[0, 0].axvline(imf.imf.mean_mass(Mmin, Mmax), **imf.plot_kwargs)

    # Cumulative number of stars
    axes[0, 1].plot(M, imf.imf.integrate(Mmin, M), label=imf.name, **imf.plot_kwargs)

    # Mass in bins
    axes[1, 0].stairs([imf.imf.integrate_product(mmin, mmax)/TOTAL_MASS for mmin, mmax in pairwise(M)], M, **imf.plot_kwargs)

    # Cumulative mass of stars
    axes[1, 1].plot(M, imf.imf.integrate_product(Mmin, M)/TOTAL_MASS, label=imf.name, **imf.plot_kwargs)

axes[0, 0].axvspan(20, 40, color="grey", alpha=0.75)
axes[0, 0].set(
    xlabel="Mass",
    ylabel="Number of stars",
    xscale="log",
    yscale="log",
    ylim=(1e-6*TOTAL_MASS, None)
    )

axes[0, 1].set(
    xlabel="Mass",
    ylabel="Cumulative Number of stars",
    xscale="log",
    yscale="log",
    ylim=(1e-4*TOTAL_MASS, None)
    )

axes[1, 0].set(
    xlabel="Mass",
    ylabel="Mass in stars/Total",
    xscale="log",
    yscale="log",
    ylim=(1e-5, None)
    )

axes[1, 1].set(
    xlabel="Mass",
    ylabel="Cumulative Mass of stars/Total",
    xscale="log",
    yscale="log",
    ylim=(1e-5, None)
    )

plt.legend()
plt.show()