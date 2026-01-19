import numpy as np
try:
    from h5py import string_dtype, File as h5File
    string_dtype = string_dtype(encoding='utf-8', length=None)
except ImportError:
    def h5File(*args, **kwargs):
        raise ImportError("The `h5py` package is required for saving to/loading from hdf5 files.")

rng = np.random.default_rng()

class IMFSample:
    """
    A class to calculate and store properties of a sample, beyond just a list of masses.

    Has the ability to store:
    * Quantity when the IMF has been sampled
    * The same Quantity when the IMF has been averaged over
    * The residual, defined as: (Qsampled - Qaveraged) / Qaveraged.
    for quantities like:
    * Number of Supernovae
    * Number of black holes
    * arbitrary quantities defined on a grid to be interpolated in mass, including ones that vary with time.

    If you want to sample from the IMF several times, calculate the same quantities, and compare them, see also IMFSampleList.

    Notes
    -----
    The residual is defined as (Qsampled - Qaveraged) / Qaveraged, so a residual of 1 means the sampled value is twice that of the averaged, or 100%.
    A value close to 0 means they are similar, far away from 0 means different.
    """
    def __init__(self, masses, target_mass, stop_method):
        self.masses = masses
        self.target_mass = target_mass
        self.total_mass = np.sum(self.masses)
        self.stop_method = stop_method
        self.sampled_quantities = {}
        self.averaged_quantities = {}
        self.residuals = {}

    def __str__(self):
        return f"Sample of {len(self.masses)} stars summing to {self.total_mass}. Sampled using `{self.stop_method}` stop_method with a target mass {self.target_mass}."

    def add_number_of_supernovae(self, imf, mass_cutoff=8):
        def count_supernovae(masses):
            return np.count_nonzero(masses > mass_cutoff)
        imf_average = imf.integrate(mass_cutoff, imf.Mmax)
        name = "Number of supernovae"
        return self.add_quantity(count_supernovae, imf_average, name)

    def add_number_of_black_holes(self, imf, mass_cutoff_lower=25, mass_cutoff_higher=100):
        def count_black_holes(masses):
            return np.count_nonzero((masses >= mass_cutoff_lower) & (masses <= mass_cutoff_higher))
        imf_average = imf.integrate(mass_cutoff_lower, mass_cutoff_higher)
        name = "Number of black holes"
        return self.add_quantity(count_black_holes, imf_average, name)

    def add_interpolated_quantity(self, mass_grid, quantity_grid, imf, name, interp_kwargs={}):
        # Provide a grid of masses and quantities, then kwargs for np.interp
        def sum_interpolated(masses):
            return np.interp(masses, mass_grid, quantity_grid, **interp_kwargs).sum()
        imf_average = imf.integrate_linear_piecewise_interpolated_product(imf.Mmin, imf.Mmax, mass_grid, quantity_grid)
        return self.add_quantity(sum_interpolated, imf_average, name)

    def add_interpolated_quantity_time_dependant(self, mass_grid, quantity_grid, imf, name, interp_kwargs={}, imf_extrapolate=False, Mmin=None, Mmax=None):
        Mmin = imf.Mmin if Mmin is None else max(Mmin, imf.Mmin)
        Mmax = imf.Mmax if Mmax is None else min(Mmax, imf.Mmax)
        # Provide a grid of masses and quantities, then kwargs for np.interp
        def sum_interpolated(masses):
            mask = (masses >= Mmin) & (masses <= Mmax)
            return np.array([np.interp(masses[mask], mass_grid, row, **interp_kwargs).sum() for row in quantity_grid])
        # I've not fully testing this broadcasting but it looks like if you do quantity_grid.T it's okay
        imf_average = imf.integrate_linear_piecewise_interpolated_product(Mmin, Mmax, mass_grid, quantity_grid.T, extrapolate_grid=imf_extrapolate)
        return self.add_quantity(sum_interpolated, imf_average, name)

    def add_quantity(self, function, imf_averaged_quantity, name):
        self.sampled_quantities[name] = function(self.masses)
        self.averaged_quantities[name] = imf_averaged_quantity
        self.residuals[name] = (self.sampled_quantities[name] - imf_averaged_quantity) / imf_averaged_quantity
        return self.sampled_quantities[name]
        # Could add getattr to read the sampled quantites dict

    def save(self, filename):
        pass  # Save to HDF5? Would be nice but adds dependencies

class IMFSampleList:
    """
    A class to calculate and store properties of several samples, beyond just a list of masses.

    Has the ability to store:
    * Quantity when the IMF has been sampled
    * The same Quantity when the IMF has been averaged over
    * The residual, defined as: (Qsampled - Qaveraged) / Qaveraged.
    for quantities like:
    * Number of Supernovae
    * Number of black holes
    * arbitrary quantities defined on a grid to be interpolated in mass, including ones that vary with time.

    If you want to sample from the IMF only once, see also IMFSample.

    Notes
    -----
    The residual is defined as (Qsampled - Qaveraged) / Qaveraged, so a residual of 1 means the sampled value is twice that of the averaged, or 100%.
    A value close to 0 means they are similar, far away from 0 means different.
    """
    def __init__(self, sample_list):
        self.sample_list = sample_list
        self.sampled_quantities = {}
        self.averaged_quantities = {}
        self.residuals = {}
        self.Nsamples = len(sample_list)
        # Need to make a decision about having different stop_methods - should it be allowed?

    def add_total_mass(self, imf):
        imf_average = imf.integrate_product(imf.Mmin, imf.Mmax)
        name = "Total Mass"
        self.add_quantity(np.sum, imf_average, name)

    def add_most_massive_star(self, imf_max=100):
        name = "Most Massive Star"
        self.add_quantity(np.max, imf_max, name)

    def add_number_of_supernovae(self, imf, mass_cutoff=8):
        def count_supernovae(masses):
            return np.count_nonzero(masses > mass_cutoff)
        imf_average = imf.integrate(mass_cutoff, imf.Mmax)
        name = "Number of supernovae"
        self.add_quantity(count_supernovae, imf_average, name)

    def add_number_of_black_holes(self, imf, mass_cutoff_lower=25, mass_cutoff_higher=100):
        def count_black_holes(masses):
            return np.count_nonzero((masses >= mass_cutoff_lower) & (masses <= mass_cutoff_higher))
        imf_average = imf.integrate(mass_cutoff_lower, mass_cutoff_higher)
        name = "Number of black holes"
        self.add_quantity(count_black_holes, imf_average, name)

    def add_interpolated_quantity(self, mass_grid, quantity_grid, imf, name, interp_kwargs={}, imf_extrapolate=False, Mmin=None, Mmax=None):
        Mmin = imf.Mmin if Mmin is None else max(Mmin, imf.Mmin)
        Mmax = imf.Mmax if Mmax is None else min(Mmax, imf.Mmax)
        # Provide a grid of masses and quantities, then kwargs for np.interp
        def sum_interpolated(masses):
            mask = (masses >= Mmin) & (masses <= Mmax)
            return np.interp(masses[mask], mass_grid, quantity_grid, **interp_kwargs).sum()
        imf_average = imf.integrate_linear_piecewise_interpolated_product(Mmin, Mmax, mass_grid, quantity_grid, extrapolate_grid=imf_extrapolate)
        return self.add_quantity(sum_interpolated, imf_average, name)

    def add_interpolated_quantity_time_dependant(self, mass_grid, quantity_grid, imf, name, interp_kwargs={}, imf_extrapolate=False, Mmin=None, Mmax=None):
        Mmin = imf.Mmin if Mmin is None else max(Mmin, imf.Mmin)
        Mmax = imf.Mmax if Mmax is None else min(Mmax, imf.Mmax)
        # Provide a grid of masses and quantities, then kwargs for np.interp
        def sum_interpolated(masses):
            mask = (masses >= Mmin) & (masses <= Mmax)
            return [np.interp(masses[mask], mass_grid, row, **interp_kwargs).sum() for row in quantity_grid]
        # I've not fully testing this broadcasting but it looks like if you do quantity_grid.T it's okay
        imf_average = imf.integrate_linear_piecewise_interpolated_product(Mmin, Mmax, mass_grid, quantity_grid.T, extrapolate_grid=imf_extrapolate)
        return self.add_quantity(sum_interpolated, imf_average, name)

    def add_quantity(self, function, imf_averaged_quantity, name):
        # We can do this list comprehension and calculate the quantity for each sample, then while IMFSample.add_quantity returns the value add it to our list at the same time
        self.sampled_quantities[name] = np.array([sample.add_quantity(function, imf_averaged_quantity, name) for sample in self.sample_list]).T
        # Take transpose so this array is (Ntime, Nsamples) to calculate residuals.
        # This way means that the first axis of all averaged, sampled, and residuals are time.
        self.averaged_quantities[name] = imf_averaged_quantity
        # Want to check if imf_averaged_quantity is an array (meaning it varies with time) or a float/int (so it doesn't)
        if isinstance(imf_averaged_quantity, np.ndarray):
            self.residuals[name] = (self.sampled_quantities[name] - imf_averaged_quantity[:, None]) / imf_averaged_quantity[:, None]
        else:
            self.residuals[name] = (self.sampled_quantities[name] - imf_averaged_quantity) / imf_averaged_quantity

    def quantile(self, name, quantiles=[0.1, 0.5, 0.9], residual=False):
        if residual is True:
            return np.quantile(self.residuals[name], quantiles)
        else:
            return np.quantile(self.sampled_quantities[name], quantiles)

    def save(self, filename):
        with h5File(filename, "w") as fileout:
            fileout.create_dataset("Number_of_samples", data=self.Nsamples)

            fileout.create_group("Samples")
            for i, sample in enumerate(self.sample_list):
                sample_group = fileout["Samples"].create_group(f"sample {i}")
                sample_group.create_dataset("Masses", data=sample.masses)
                sample_group["Masses"].attrs.create("Units", "solar mass", dtype=string_dtype)
                sample_group.create_dataset("Target_mass", data=sample.target_mass)
                sample_group["Target_mass"].attrs.create("Units", "solar mass", dtype=string_dtype)
                sample_group.create_dataset("Stop_method", data=sample.stop_method, dtype=string_dtype)

            fileout.create_group("Derived_quantities")
            for quantity in set(self.averaged_quantities).union(self.sampled_quantities, self.residuals):
                quantity_group = fileout["Derived_quantities"].create_group(quantity)
                if quantity in self.averaged_quantities:
                    quantity_group.create_dataset("IMF Averaged", data=self.averaged_quantities[quantity])
                if quantity in self.sampled_quantities:
                    quantity_group.create_dataset("Sampled", data=self.sampled_quantities[quantity])
                if quantity in self.residuals:
                    quantity_group.create_dataset("Residuals", data=self.residuals[quantity])
                    quantity_group["Residuals"].attrs.create("Definition", "(sampled_quantity - averaged_quantity) / averaged_quantity")

    @classmethod
    def load(cls, filename):
        with h5File(filename, "r") as filein:
            Nsamples = filein["Number_of_samples"][()]

            sample_list = []
            for i in range(Nsamples):
                sample_group = filein[f"Samples/sample {i}"]
                masses = sample_group["Masses"][:]
                target_mass = sample_group["Target_mass"][()]
                stop_method = sample_group["Stop_method"].asstr()[()]
                sample_list.append(IMFSample(masses, target_mass, stop_method))

            samples = cls(sample_list)

            for quantity in filein["Derived_quantities"]:
                samples.averaged_quantities[quantity] = filein[f"Derived_quantities/{quantity}/IMF Averaged"][()]
                samples.sampled_quantities[quantity] = filein[f"Derived_quantities/{quantity}/Sampled"][()]
                samples.residuals[quantity] = filein[f"Derived_quantities/{quantity}/Residuals"][()]

            return samples

def draw_samples(imf, stop_method="below", target_mass=None, full_output=False, rescale=False, rng=rng):
    """
    _summary_

    Parameters
    ----------
    imf : InitialMassFunction subclass
        The IMF to sample from. Needs to have `inverse_cdf()` method implemented.
    stop_method : str from {"below" | "above" | "closest"}, optional
        The stopping criteria, see notes below for more information, by default "below"
    target_mass : int or float, optional
        The target mass to aim for, by default None. If None, use the total mass of the provided IMF.
    full_output : bool, optional
        Whether or not to return an IMFSample object (True) or just a numpy array of the sampled masses, by default False (only return mass array).
    rescale : bool, optional
        Whether or not to rescale each mass so the final total is the target mass exactly. mi -> mi * Mtarget / Mtotal. By default, False (don't rescale).
    rng : instance on numpy random Generator, optional
        The instance of numpy Generator to use. This is how the user can provide their own seed/random number generation method. By default the result of calling numpy.random.default_rng().

    Returns
    -------
    ndarray or IMFSample object
        An array of the masses drawn randomly from the IMF if full_output=False (default), or an IMFSample object representing the sample and additional information.

    """
    if stop_method not in (stop_options:=["below", "above", "closest"]):
        raise ValueError(f"`{stop_method=}` is not a valid argument. Choose from {stop_options}.")

    # Let people either declare target_mass, or we can integrate over the imf provided.
    # This means we will be able to implement Smith (2021) sampling scheme: https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.5417S/abstract
    # As here the target_mass changes.
    target_mass = imf.integrate_product(imf.Mmin, imf.Mmax) if target_mass is None else target_mass
    Nguess = 1 + 2*int(imf.integrate_product(imf.Mmin, imf.Mmax) * target_mass / imf.integrate_product(imf.Mmin, imf.Mmax))  # average number of stars, plus 1 in case someone puts something silly as an input and we end up in an infinite loop.
    drawn_masses = imf.inverse_cdf(rng.uniform(size=Nguess))
    Mcurrent = np.cumsum(drawn_masses)

    while np.all(Mcurrent < target_mass):
        # Draw more masses - Could be more intelligent about how we pick Nguess
        drawn_masses = np.append(drawn_masses, imf.inverse_cdf(rng.uniform(size=Nguess)))
        Mcurrent = np.cumsum(drawn_masses)  # Potentially can find better way to do this

    # Target hit, move onto logic
    # find the index where we are ABOVE the target.
    i = (Mcurrent > target_mass).nonzero()[0][0]  # .nonzero() will raise index error if can't find anything

    if stop_method == "below":
        sample = drawn_masses[:i]
    elif stop_method == "above":
        sample = drawn_masses[:i+1]
    elif stop_method == "closest":
        # 'Draw' the mass and add it to the list if it gets us closer to the target than leaving it off would
        # ('Draw' because we have already technically drawn it, but you know what I mean)
        if (Mcurrent[i] - target_mass) < (target_mass - Mcurrent[i-1]):
        # if drawn_mass < 2 * (target_mass - Mcurrent):
            sample = drawn_masses[:i+1]
        else:
            sample = drawn_masses[:i]

    if rescale:
        sample *= target_mass / sample.sum()  # Rescale each sampled mass so that the total mass is correct
        # Should have some check/clipping that we don't go outside of our mass range

    if full_output:
        return IMFSample(sample, target_mass, stop_method)
    else:
        return sample

def smith2021_sampling(imf, Nsamples, Mtarget=None):
    # The paper uses the terms
    # start with Moffset = 0
    # Select the particle to populate with Mpart
    # Mtarget = Mpart - Moffset
    # draw masses until Masn > Mtarget (what I would call "above" stop condition)
    # set new Moffset = Masn - Mtarget
    Mtarget = imf.integrate_product(imf.Mmin, imf.Mmax) if Mtarget is None else Mtarget
    Moffset = 0
    samples = []
    for _ in range(Nsamples):
        sample = draw_samples(imf, stop_method="above", target_mass=(Mtarget - Moffset))
        samples.append(sample)
        Moffset = sample.sum() - (Mtarget - Moffset)
    return samples