import os
import os.path as pth

from numpy.distutils.core import Extension

try:
    adashome = os.environ['ADASHOME']
except KeyError:
    adashome = os.path.expanduser('~adas')

os.environ['FC'] = 'gfortran'
os.environ['F77'] = 'gfortran'
os.environ['FFLAGS'] = '-fPIC -fno-automatic'
extra_link_args_f2py = ['-shared']
extra_f77_compile_args = []

files_c5dplr_ = ['fortran/adas3xx/adas305/c5dplr.for',
                 'fortran/adas3xx/adas305/c5pixv.for',
                 'fortran/adaslib/maths/r8erfc.for']
files_c5dplr_ = [pth.join(adashome, x) for x in files_c5dplr_]
files_c5dplr_.insert(0, 'cherab/adas/beam_emission/stark/c5dplr_.pyf')

files_stark_ = ['fortran/adas3xx/adas305/stark.for',
                'fortran/adas2xx/adas211/delta2.for',
                'fortran/adas3xx/adas305/bornp1.for',
                'fortran/adas3xx/adas305/bornp2.for',
                'fortran/adas3xx/adas305/c5rlsp.for',
                'fortran/adas3xx/adas305/drv.for',
                'fortran/adas3xx/adas305/drvspl.for',
                'fortran/adas3xx/adas305/fsplin.for',
                'fortran/adas3xx/adas305/hydemi.for',
                'fortran/adas3xx/adas305/lumsis.for',
                'fortran/adas3xx/adas305/sfi2.for',
                'fortran/adas3xx/adas305/sigma.for',
                'fortran/adas3xx/adas305/sigmel.for',
                'fortran/adas3xx/adas305/sigmin.for',
                'fortran/adas3xx/adas305/stark2.for',
                'fortran/adas3xx/adas305/unbun2.for',
                'fortran/adas3xx/adas305/watint.for',
                'fortran/adas3xx/adas305/zeemn2.for',
                'fortran/adas3xx/adas311/gamaf.for',
                'fortran/adas3xx/adas311/psprod.for',
                'fortran/adas3xx/adas311/rwfh.for',
                'fortran/adas8xx/adas804/dipol.for',
                'fortran/adaslib/atomic/wig3j.for',
                'fortran/adaslib/system/i4unit.for',
                'fortran/adaslib/maths/xxminv.for']

files_stark_ = [pth.join(adashome, x) for x in files_stark_]
files_stark_.insert(0, 'cherab/adas/beam_emission/stark/stark_.pyf')

ext_modules = []

libs_stark_ = ['lapack', 'blas']
extra_link_args_stark_ = list(extra_link_args_f2py)  # copy


c5dplr_ = Extension(name='cherab.adas.beam_emission.stark.c5dplr_',
                    sources=files_c5dplr_,
                    extra_f77_compile_args=extra_f77_compile_args,
                    extra_link_args=extra_link_args_f2py)

stark_ = Extension(name='cherab.adas.beam_emission.stark.stark_',
                   sources=files_stark_,
                   extra_f77_compile_args=extra_f77_compile_args,
                   extra_link_args=extra_link_args_stark_,
                   libraries=libs_stark_)

ext_modules += [c5dplr_, stark_]

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(ext_modules=ext_modules)
