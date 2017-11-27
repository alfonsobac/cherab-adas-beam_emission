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

from raysect.optical cimport Spectrum, Point3D, Vector3D

from cherab.core cimport Line, AtomicData
from cherab.core.beam cimport BeamModel


cdef class BeamEmissionLine(BeamModel):

    cdef:
        Line _line
        double _wavelength
        double _plasma_mass
        double _beam_temperature

        double _pi_factor
        double _sigma_factor
        list _beam_emission_data

    cdef inline Spectrum _beam_line_emission(self, Spectrum spectrum, double x, double y, double z,
                                             double beam_density, Vector3D observation_direction, Vector3D beam_direction)

    cdef inline double _beam_emission(self, double x, double y, double z, Vector3D beam_velocity)

    cdef inline int _populate_cache(self) except -1
