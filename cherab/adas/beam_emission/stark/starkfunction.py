import numpy as np

from cherab.core.utility import PerCm3ToPerM3, AngstromToNm

from .c5dplr_ import c5dplr as _c5dplr_fort
from .stark_ import stark as _stark_fort

STARK_D_ALPHA_WAVELENGTH = 656.09522


def stark(beam_mass, beam_energy, beam_temperature,
          beam_density, beam_direction,
          plasma_mass, plasma_electron_temperature,
          plasma_electron_density, plasma_zeff,
          b_field, observation_direction, sigma_factor, pi_factor, transition,
          wavelength_size, wavelength_min, wavelength_max,
          d_alpha_wavelength=0.0):
    """Calculates the components of the Hydrogen Stark feature with
    its Doppler broadening.

    :param beam_mass: H beam mass (amu)
    :param beam_energy: In eV/amu
    :param beam_temperature: In eV
    :param beam_density: In :math:`m^{-3}`
    :param beam_direction: Vector :math:`(\\hat{v}_x, \\hat{v}_y, \\hat{v}_z)`
    :param plasma_mass: H mass in plasma (amu)
    :param plasma_electron_temperature: In eV
    :param plasma_electron_density: In :math:`m^{-3}`
    :param plasma_zeff: Effective charge :math:`Z_{eff}`
    :param b_field: magnetic field vector (3-component vector in Tesla)
    :param observation_direction:  Normalized vector
    :param sigma_factor: :math:`\\sigma` polarisation transmission factor
    :param pi_factor: :math:`\\pi` polarisation transmission factor
    :param transition: Tuple with: (initial_level, final_level)
    :param wavelength_size: actual output array number of elements
    :param wavelength_min: minimum wavelength for output array (nm)
    :param wavelength_max: maximum wavelength for output array (nm)
    :returns: array of floats -- Normalised Doppler broadened feature.
              The sum of returned array is 1
    :param d_alpha_wavelength: Real D_alpha wavelength; if <= 0.0,
                               there is not correction in wavelenght
    """

    if d_alpha_wavelength > 1.0:
        wavelength_shift = d_alpha_wavelength - STARK_D_ALPHA_WAVELENGTH
    else:
        wavelength_shift = 0.0

    # From nm to angstrom
    wvmin_a = AngstromToNm.inv(wavelength_min - wavelength_shift)
    wvmax_a = AngstromToNm.inv(wavelength_max - wavelength_shift)

    # Wavelegths will be angstrom from here

    amdeut = beam_mass
    bener = beam_energy
    tbev = beam_temperature
    dv1, dv2, dv3 = beam_direction
    densb = PerCm3ToPerM3.inv(beam_density)   # Convert from m-3 to cm-3

    amss = plasma_mass
    te = plasma_electron_temperature
    zeff = plasma_zeff
    dens = PerCm3ToPerM3.inv(plasma_electron_density)  # From m-3 to cm-3

    bmag = b_field.length
    if bmag > 0:
        db1, db2, db3 = b_field / bmag
    else:
        db1, db2, db3 = 1.0, 0.0, 0.0

    de1, de2, de3 = 1.0, 0.0, 0.0
    emag = 0.0                             # V/cm

    # Change direction to be right
    do1, do2, do3 = - observation_direction
    polo = sigma_factor
    polp = pi_factor

    nu, nl = transition

    popu = 1.0

    # set to 300 the max of components (not more)
    ndcomp = 300

    wvcomp = np.zeros(ndcomp)
    emcomp = np.zeros(ndcomp)

    polofull = 1.0
    polpfull = 1.0

    # Calculate stark with full transmission

    ncomp, _, _ = _stark_fort(amdeut, amss, bener, dv1, dv2, dv3, densb,
                              bmag, db1, db2, db3, emag, de1, de2, de3,
                              do1, do2, do3, polofull, polpfull,
                              dens, te, zeff, nu, nl, popu, wvcomp, emcomp)

    doppler = np.zeros(wavelength_size)

    sumcomp = emcomp.sum()  # Normalization

    if ncomp == 0 or sumcomp <= 0:
        return doppler

    # If transmission is not full, it should be calculated again
    if polo != 1.0 or polp != 1.0:
        wvcomp = np.zeros(ndcomp)
        emcomp = np.zeros(ndcomp)

        ncomp, _, _ = _stark_fort(amdeut, amss, bener, dv1, dv2, dv3, densb,
                                  bmag, db1, db2, db3, emag, de1, de2, de3,
                                  do1, do2, do3, polo, polp,
                                  dens, te, zeff, nu, nl, popu, wvcomp, emcomp)

        if ncomp == 0:
            return doppler

    wvcomp = wvcomp[:ncomp]
    emcomp = emcomp[:ncomp]

    w_min = wvcomp.min()
    w_max = wvcomp.max()

    extra_band = 1.4e-4 * w_max * np.sqrt(tbev / amdeut)
    w_min -= extra_band
    w_max += extra_band

    if wvmax_a < w_min or wvmin_a > w_max:
        # If wvmin and wvmax are far from w_min and w_max,
        # a core dump can be produced
        return doppler

    _c5dplr_fort(wavelength_size, wvmin_a, wvmax_a, ncomp, wvcomp, emcomp,
                 tbev, amdeut, doppler)

    # doppler.sum() == emcomp.sum(), but emcomp always includes full range

    return doppler / sumcomp
