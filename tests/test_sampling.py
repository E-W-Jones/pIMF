import numpy as np
import matplotlib.pyplot as plt
import pimf

implemented_imfs = {
    "Power Law": pimf.PowerLawIMF(normalisation="mass", normalisation_value=1e5),
    "Broken Power Law": pimf.BrokenPowerLawIMF(normalisation="mass", normalisation_value=1e5),
    "Chabrier": pimf.ChabrierIMF(normalisation="mass", normalisation_value=1e5),
    "Lognormal": pimf.LognormalIMF(normalisation="mass", normalisation_value=1e5, mc=10, Mmin=0.1, Mmax=100),
    "L3": pimf.L3IMF(normalisation="mass", normalisation_value=1e5),
    "Generalised Gamma": pimf.GeneralisedGammaIMF(normalisation="mass", normalisation_value=1e6, Mmin=0.1, Mmax=100, mc=1)
}

def test_inverse_CDF():
    # The idea of this is to put in a mass, get a CDF for that mass, then plug it back through the inverse to get the same mass.
    pass

def test_imf_visually():
    M = np.geomspace(0.1, 100)
    fig, axes = plt.subplots(3, 2, sharex=True, sharey=True)
    for ax, imf_name in zip(axes.flatten(), implemented_imfs):
        imf = implemented_imfs[imf_name]
        ax.loglog(M, imf(M) / imf.integrate(0.1, 100), label="IMF")
        ax.hist(pimf.draw_samples(imf), bins=M, density=True, label="Drawn Masses")
        ax.set_title(imf_name)
        ax.set_ylim(9e-7)
    fig.supxlabel("Mass / Msol")
    fig.supylabel("PDF")
    axes[0, 0].legend()
    plt.show()