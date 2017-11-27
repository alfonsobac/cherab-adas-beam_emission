# cython: language_level=3

# Copyright 2014-2017 United Kingdom Atomic Energy Authority
#
# Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the
# European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/software/page/eupl5
#
# Unless required by applicable law or agreed to in writing, software distributed
# under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied.
#
# See the Licence for the specific language governing permissions and limitations
# under the Licence.

"""Calculate Beam emission with ADAS beam rates and stark routine"""

from scipy import constants

from libc.math cimport exp, sqrt, M_PI as pi
from numpy cimport ndarray
cimport cython
from raysect.optical.material.emitter.inhomogeneous import NumericalIntegrator
from raysect.optical cimport Spectrum, Point3D, Vector3D

from cherab.core cimport Species, Plasma, Beam, Element, BeamEmissionRate
from cherab.core.model.lineshape import doppler_shift, thermal_broadening, add_gaussian_line
from cherab.core.utility.constants cimport RECIP_4_PI, ELEMENTARY_CHARGE, ATOMIC_MASS
from cherab.adas.beam_emission.stark import stark

cdef double RECIP_ELEMENTARY_CHARGE = 1 / ELEMENTARY_CHARGE
cdef double RECIP_ATOMIC_MASS = 1 / ATOMIC_MASS


cdef double evamu_to_ms(double x):
    return sqrt(2 * x * ELEMENTARY_CHARGE * RECIP_ATOMIC_MASS)


cdef double ms_to_evamu(double x):
    return 0.5 * (x ** 2) * RECIP_ELEMENTARY_CHARGE * ATOMIC_MASS


cdef class BeamEmissionLine(BeamModel):
    """Calculates Beam emission from a beam (collisional and Stark).

    :param line: 
    :param beam: 
    :param plasma: 
    :param atomic_data: 
    :param sigma_factor: factor for :math:`\\sigma` polarised light
    :param pi_factor: factor for :math:`\\pi` polarised light
    :param plasma_mass: H plasma mass for ADAS stark function
    :param beam_temperature: beam temperature (eV) for ADAS stark function   
    :return:
    """

    def __init__(self, Line line not None, 
                 Beam beam=None, Plasma plasma=None, AtomicData atomic_data=None,
                 double sigma_factor=1.0, double pi_factor=1.0, 
                 double plasma_mass=0, double beam_temperature=0):

        super().__init__(beam, plasma, atomic_data)

        self.line = line
        self.plasma_mass = plasma_mass
        self.beam_temperature = beam_temperature
        self.sigma_factor = sigma_factor
        self.pi_factor = pi_factor


        # initialise cache to empty
        self._change()


    @property
    def line(self):
        return self._line

    @line.setter
    def line(self, Line value not None):
        # the data cache depends on the line configuration
        self._line = value
        self._change()

    @property    
    def plasma_mass(self):
        return self._plasma_mass

    @plasma_mass.setter
    def plasma_mass(self, double plasma_mass):
        if plasma_mass <= 0:
            raise ValueError("Plasma mass has to be bigger than 0")
        self._plasma_mass = plasma_mass

    @property    
    def beam_temperature(self):
        return self._beam_temperature

    @beam_temperature.setter
    def beam_temperature(self, double beam_temperature):
        if beam_temperature <= 0:
            raise ValueError("Beam temperature has to be bigger than 0")
        self._beam_temperature = beam_temperature

    @property    
    def pi_factor(self):
        return self._pi_factor

    @pi_factor.setter
    def pi_factor(self, double pi_factor):
        if pi_factor < 0:
            raise ValueError("Pi factor can not be less than 0")
        self._pi_factor = pi_factor

    @property    
    def sigma_factor(self):
        return self._sigma_factor

    @sigma_factor.setter
    def sigma_factor(self, double sigma_factor):
        if sigma_factor < 0:
            raise ValueError("Pi factor can not be less than 0")
        self._sigma_factor = sigma_factor


    # todo: escape early if data is not suitable for a calculation
    # todo: carefully review changes to maths
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cpdef Spectrum emission(self, Point3D beam_point, Point3D plasma_point, Vector3D beam_direction,
                            Vector3D observation_direction, Spectrum spectrum):

        cdef:
            double beam_density, x, y, z

        # cache data on first run
        if self._beam_emission_data is None:
            self._populate_cache()

        # obtain beam density
        beam_density = self._beam.density(beam_point.x, beam_point.y, beam_point.z)

        # abort calculation if beam density is zero
        if beam_density == 0.0:
            return spectrum


        # extract for more compact code
        x = plasma_point.x
        y = plasma_point.y
        z = plasma_point.z


        return self._beam_line_emission(spectrum, x, y, z, beam_density, observation_direction, beam_direction)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline Spectrum _beam_line_emission(self, Spectrum spectrum, double x, double y, double z,
                                            double beam_density, Vector3D observation_direction, Vector3D beam_direction):
        """
        Beam emission spectrum calculated with ADAS beam rates

        :param spectrum: Raysect spectrum object
        :param x: position in meters
        :param y: position in meters
        :param z: position in meters
        :param beam_density: neutral beam density in m^-3
        :param observation_direction:
        :param beam_direction:
        :return: spectral line emission in W/m^3/str/nm
        """
        cdef:
            double beam_emission, delta_wavelength, wavelength_min, wavelength_max
            double beam_energy, beam_mass, beam_temperature
            double plasma_zeff, plasma_electron_density, plasma_electron_temperature
            double plasma_mass
            Vector3D beam_velocity, b_field
            int wavelength_size
            tuple transition
            double sigma_factor, pi_factor

        beam_velocity = beam_direction * evamu_to_ms(self._beam.energy)
        beam_emission = self._beam_emission(x, y, z, beam_velocity)
        
                
        wavelength_size = spectrum.bins
        delta_wavelength = spectrum.delta_wavelength
        
        wavelength_min = spectrum.wavelengths[0]
        wavelength_max = spectrum.wavelengths[wavelength_size - 1]

        beam_energy = self._beam.energy
        beam_mass = self._beam.element.atomic_weight
        beam_temperature = self.beam_temperature
        
        # calculate z_effective and the B-field magnitude
        try:
            plasma_zeff = self._plasma.z_effective(x, y, z)
        except ValueError:
            # There is no plasma or no species
            return spectrum

        b_field = self._plasma.get_b_field().evaluate(x, y, z)
        
        plasma_electron_density = self._plasma.get_electron_distribution().density(x, y, z)
        plasma_electron_temperature = self._plasma.get_electron_distribution().effective_temperature(x, y, z)
        
        plasma_mass = self.plasma_mass
        transition = self._line.transition
        sigma_factor = self.sigma_factor
        pi_factor = self.pi_factor
        
        stark_features = stark(beam_mass, beam_energy, beam_temperature, 
                               beam_density, beam_direction, 
                               plasma_mass, plasma_electron_temperature, 
                               plasma_electron_density, plasma_zeff,
                               b_field, observation_direction, sigma_factor, pi_factor, 
                               transition, wavelength_size, wavelength_min, wavelength_max,
                               d_alpha_wavelength=self._wavelength)
        
        spectrum.samples[:] += stark_features * beam_emission * beam_density / (4 * pi * delta_wavelength)

        return spectrum



    @cython.cdivision(True)
    cdef inline double _beam_emission(self, double x, double y, double z, Vector3D beam_velocity):
        """
        
        :param x: position in meters
            return spectrum
        :param y: position in meters
        :param z: position in meters


        :param beam_velocity: beam velocity in m/s
        :return: a beam emission rate in s^-1
        """

        # see www.adas.ac.uk/man/chap3-04.pdf equation 4.4.7
        # note: we have access to ni for each species so we have done away with
        # the impurity fractions used in the above document

        cdef:
            double density_sum, beam_emission, target_ne, target_ti, interaction_speed, interaction_energy, target_equiv_ne
            Species species
            BeamEmissionRate be_rate
            int target_z
            Vector3D target_velocity, interaction_velocity

        # z-weighted density sum
        density_sum = 0
        for species, _ in self._beam_emission_data:
            density_sum += species.element.atomic_number**2 * species.distribution.density(x, y, z)

        # beam emission rate
        beam_emission = 0
        for species, be_rate in self._beam_emission_data:

            # sample species distribution
            target_z = species.element.atomic_number
            target_ne = species.distribution.density(x, y, z) * target_z
            target_ti = species.distribution.effective_temperature(x, y, z)
            target_velocity = species.distribution.bulk_velocity(x, y, z)

            # calculate mean beam interaction energy
            interaction_velocity = beam_velocity - target_velocity
            interaction_speed = interaction_velocity.length
            interaction_energy = ms_to_evamu(interaction_speed)

            # species equivalent electron density
            target_equiv_ne = density_sum / target_z

            beam_emission += target_ne * be_rate.evaluate(interaction_energy, target_equiv_ne, target_ti)

        return beam_emission

    cdef inline int _populate_cache(self) except -1:

        cdef:
            Element target_element, beam_element
            tuple transition
            Species species
            BeamEmissionRate be_rate

        # sanity checks
        if self._beam is None or self._plasma is None or self._atomic_data is None:
            raise RuntimeError("The emission model is not connected to a beam object.")

        if self._line is None:
            raise RuntimeError("The emission line has not been set.")

        target_element = self._line.element
        beam_element = self._beam.element

        if target_element != beam_element:
            raise RuntimeError("Target element should be beam element.")

        transition = self._line.transition

        # obtain wavelength for specified line
        self._wavelength = self._atomic_data.wavelength(target_element, 0, transition)

        # obtain beam emission rates
        self._beam_emission_data = []
        for species in self._plasma.composition:
            be_rate = self._atomic_data.beam_emission_rate(beam_element, species.element, species.ionisation, transition)
            self._beam_emission_data.append((species, be_rate))
        

    def _change(self):
        # clear cache to force regeneration on first use
        self._beam_emission_data = None
        self._wavelength = 0.0
