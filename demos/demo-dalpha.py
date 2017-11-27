from scipy.constants import electron_mass, atomic_mass

from cherab.core import Plasma, Beam, Species, Maxwellian, elements, Line
from cherab.openadas import OpenADAS
from cherab.core.model import BeamCXLine, SingleRayAttenuator
from cherab.adas.beam_emission import BeamEmissionLine
from raysect.optical import World, Point3D, Vector3D, Ray, translate, rotate
from cherab.core.math import ConstantVector3D
from gaussian_volume import GaussianVolume

integration_step = 0.02

# setup scenegraph
world = World()

# create atomic data source
adas = OpenADAS(permit_extrapolation=True)

# PLASMA ----------------------------------------------------------------------
plasma = Plasma(parent=world)


# define basic distributions
ion_density = 5e19
sigma = 0.25


d_density = GaussianVolume(0.94 * ion_density, sigma)
he2_density = GaussianVolume(0.04 * ion_density, sigma)
c6_density = GaussianVolume(0.01 * ion_density, sigma)
ne10_density = GaussianVolume(0.01 * ion_density, sigma)
e_density = GaussianVolume((0.94 + 0.04*2 + 0.01*6 + 0.01*10) * ion_density, sigma)
temperature = 1000 + GaussianVolume(4000, sigma)
bulk_velocity = ConstantVector3D(Vector3D(200e3, 0, 0))


d_distribution = Maxwellian(d_density, temperature, bulk_velocity, 
                            elements.deuterium.atomic_weight * atomic_mass)
he2_distribution = Maxwellian(he2_density, temperature, bulk_velocity, 
                              elements.helium.atomic_weight * atomic_mass)
c6_distribution = Maxwellian(c6_density, temperature, bulk_velocity, 
                             elements.carbon.atomic_weight * atomic_mass)
ne10_distribution = Maxwellian(ne10_density, temperature, bulk_velocity, 
                               elements.neon.atomic_weight * atomic_mass)
e_distribution = Maxwellian(e_density, temperature, bulk_velocity, 
                            electron_mass)

d_species = Species(elements.deuterium, 1, d_distribution)
he2_species = Species(elements.helium, 2, he2_distribution)
c6_species = Species(elements.carbon, 6, c6_distribution)
ne10_species = Species(elements.neon, 10, ne10_distribution)

# define speciess.samples.max()
plasma.b_field = ConstantVector3D(Vector3D(0.0, 0.0, 3.0))
plasma.electron_distribution = e_distribution
plasma.composition = [d_species, he2_species, c6_species, ne10_species]

# EMISSION MODEL

deuteriumline = Line(elements.deuterium, 0, (3, 2))
deuterium_emission = BeamCXLine(deuteriumline)

be_emission = BeamEmissionLine(deuteriumline,
                                plasma_mass=elements.deuterium.atomic_weight,
                                sigma_factor=1.0, pi_factor=1.0,
                                beam_temperature=5)
emission_models = [deuterium_emission, be_emission]


# BEAM ------------------------------------------------------------------------
beam = Beam(parent=world, transform=translate(2.0, 0.0, 0.5) * rotate(120, 0, 0))
beam.plasma = plasma
beam.atomic_data = adas
beam.energy = 40000
beam.power = 3e6
beam.element = elements.deuterium
beam.sigma = 0.025
beam.divergence_x = 0.5
beam.divergence_y = 0.5
beam.length = 3.0
beam.attenuator = SingleRayAttenuator(clamp_to_zero=True)
beam.models = emission_models
beam.integrator.step = integration_step
beam.integrator.min_samples = 10


# OBSERVER --------------------------------------------------------------------

# Ray goes to z-axis (positive)
r = Ray(origin=Point3D(0, 0, -3), min_wavelength=650, max_wavelength=660,
        bins=1000)

s = r.trace(world)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(s.wavelengths, s.samples)
plt.show()
