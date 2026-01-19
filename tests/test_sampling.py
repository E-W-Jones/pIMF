import numpy as np
import matplotlib.pyplot as plt
import pimf

implemented_imfs = {
    "Power Law": pimf.PowerLawIMF(normalisation="mass", normalisation_value=1e5),
    "Power Law (alpha=-1)": pimf.PowerLawIMF(alpha=-1, normalisation="mass", normalisation_value=1e5),
    "Broken Power Law": pimf.BrokenPowerLawIMF(normalisation="mass", normalisation_value=1e5),
    "Chabrier": pimf.ChabrierIMF(normalisation="mass", normalisation_value=1e5),
    "Lognormal": pimf.LognormalIMF(normalisation="mass", normalisation_value=1e5, mc=10, Mmin=0.1, Mmax=100),
    "L3": pimf.L3IMF(normalisation="mass", normalisation_value=1e5),
    "Generalised Gamma": pimf.GeneralisedGammaIMF(normalisation="mass", normalisation_value=1e6, Mmin=0.1, Mmax=100, mc=1)
}

implemented_imfs_low_mass = {
    "Power Law": pimf.PowerLawIMF(normalisation="mass", normalisation_value=1e2),
    "Power Law (alpha=-1)": pimf.PowerLawIMF(alpha=-1, normalisation="mass", normalisation_value=1e2),
    "Broken Power Law": pimf.BrokenPowerLawIMF(normalisation="mass", normalisation_value=1e2),
    "Chabrier": pimf.ChabrierIMF(normalisation="mass", normalisation_value=1e2),
    "Lognormal": pimf.LognormalIMF(normalisation="mass", normalisation_value=1e2, mc=10, Mmin=0.1, Mmax=100),
    "L3": pimf.L3IMF(normalisation="mass", normalisation_value=1e2),
    "Generalised Gamma": pimf.GeneralisedGammaIMF(normalisation="mass", normalisation_value=1e2, Mmin=0.1, Mmax=100, mc=1)
}

def test_inverse_CDF():
    # The idea of this is to put in a mass, get a CDF for that mass, then plug it back through the inverse to get the same mass.
    pass

def test_imf_visually():
    M = np.geomspace(0.1, 100)
    fig, axes = plt.subplots(4, 2, sharex=True, sharey=True)
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

def test_interpolation():
    # pick a really simple grid, just using L propto M^3 in solar units
    mass_grid = np.geomspace(0.1, 100, 10)
    L_grid = mass_grid**3
    for name, imf in implemented_imfs_low_mass.items():
        sample = pimf.draw_samples(imf, full_output=True)
        print(sample)
        sample.add_interpolated_quantity(
            mass_grid=mass_grid, quantity_grid=L_grid, imf=imf, name="Luminosity"
        )
        print(sample.averaged_quantities["Luminosity"], sample.sampled_quantities["Luminosity"])

def test_interpolation_over_time_manual():
    # pick a really simple grid, just using L propto M^3 in solar units
    # But this time we will also have a star only emit for some if its life, following t propto M^-2 (L propto M/t propto M^3)
    NMass = 10
    Ntime = 20
    mass_grid = np.geomspace(0.1, 100, NMass)
    t_grid = np.geomspace(110**-2, 0.09**-2, Ntime)
    L_grid = np.tile(mass_grid**3, (Ntime, 1)) * (t_grid[:, None] < mass_grid[None, :]**-2)

    fig, axes = plt.subplots(4, 2, sharex=True, sharey=True)

    for ax, imf_name in zip(axes.flatten(), implemented_imfs_low_mass):
        imf = implemented_imfs[imf_name]
        sample = pimf.draw_samples(imf, full_output=True)
        # print(sample)
        for i in range(Ntime):
            sample.add_interpolated_quantity(
                mass_grid=mass_grid, quantity_grid=L_grid[i], imf=imf, name=f"Luminosity {t_grid[i]}"
            )
            # print(sample.averaged_quantities[f"Luminosity {t_grid[i]}"], sample.sampled_quantities[f"Luminosity {t_grid[i]}"])
        ax.loglog(t_grid, [sample.sampled_quantities[f"Luminosity {t}"] for t in t_grid], label="Interpolated")
        ax.plot(t_grid, 1e4*t_grid**-0.825, label="Analytical $L\\propto t^{-0.825}$")
        ax.set_title(imf_name)
    fig.supxlabel("time / solar lifetimes")
    fig.supylabel("Luminosity / solar units")
    axes[0, 0].legend()
    plt.show()

def test_interpolation_over_time_function():
    # pick a really simple grid, just using L propto M^3 in solar units
    # But this time we will also have a star only emit for some if its life, following t propto M^-2 (L propto M/t propto M^3)
    NMass = 10
    Ntime = 20
    mass_grid = np.geomspace(0.1, 100, NMass)
    t_grid = np.geomspace(110**-2, 0.09**-2, Ntime)
    L_grid = np.tile(mass_grid**3, (Ntime, 1)) * (t_grid[:, None] < mass_grid[None, :]**-2)

    fig, axes = plt.subplots(4, 2, sharex=True, sharey=True)

    for ax, imf_name in zip(axes.flatten(), implemented_imfs_low_mass):
        imf = implemented_imfs[imf_name]
        sample = pimf.draw_samples(imf, full_output=True)
        # print(sample)
        sample.add_interpolated_quantity_time_dependant(
            mass_grid=mass_grid, quantity_grid=L_grid, imf=imf, name="Luminosity"
        )
        ax.loglog(t_grid, sample.sampled_quantities["Luminosity"])
        ax.set_title(imf_name)
    fig.supxlabel("time / solar lifetimes")
    fig.supylabel("Luminosity / solar units")
    # axes[0, 0].legend()
    plt.show()

def test_interpolation_over_time_function_samples():
    # pick a really simple grid, just using L propto M^3 in solar units
    # But this time we will also have a star only emit for some if its life, following t propto M^-2 (L propto M/t propto M^3)
    NMass = 10
    Ntime = 20
    mass_grid = np.geomspace(0.1, 100, NMass)
    t_grid = np.geomspace(110**-2, 0.09**-2, Ntime)
    L_grid = np.tile(mass_grid**3, (Ntime, 1)) * (t_grid[:, None] < mass_grid[None, :]**-2)

    fig, axes = plt.subplots(4, 2, sharex=True, sharey=True)

    for ax, imf_name in zip(axes.flatten(), implemented_imfs_low_mass):
        imf = implemented_imfs[imf_name]
        sample = pimf.sampling.IMFSampleList([pimf.draw_samples(imf, full_output=True) for _ in range(2)])
        # print(sample)
        sample.add_interpolated_quantity_time_dependant(
            mass_grid=mass_grid, quantity_grid=L_grid, imf=imf, name="Luminosity"
        )
        # Plot the residual as we have a large dynamic range
        ax.semilogx(t_grid, sample.residuals["Luminosity"], color="grey")
        ax.set_title(imf_name)
    fig.supxlabel("time / solar lifetimes")
    fig.supylabel("Luminosity / solar units residual")
    # axes[0, 0].legend()
    plt.show()