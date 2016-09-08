from distutils.core import setup

setup(
    name='t80s_reduce',
    version='0.0.1',
    packages=['t80s_reduce', 'fits_reduce.main', 'fits_reduce.util'],
    url='https://github.com/tribeiro/t80s_reduce/',
    license='BSD',
    author='Tiago Ribeiro de Souza',
    author_email='tribeiro@ufs.br',
    description='An image reducer for robotic telescopes',
    install_requires=['astropy',
                      'colorlog',
                      'ccdproc>=0.3.3',
                      'fit_reduce'],
    scripts=['scripts/reducer',
             'scripts/t80s_imarith',
             'scripts/t80s_imcombine',
             'scripts/t80s_preproc',]
)
