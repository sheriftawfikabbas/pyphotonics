from distutils.core import setup
setup(
    name='pyphotonics',
    packages=['pyphotonics'],
    version='0.1.6',
    license='gpl-3.0',
    description='The PyPhotonics python code is a post-processing code written entirely in python which takes as input the output files of the VASP and phonopy codes for a defect system, and calculates the Huang-Rhys factor and the PL lineshapes for that system.',   # Give a short description about your library
    author='Sherif Abdulkader Tawfik',
    author_email='sherif.tawfic@gmail.com',
    url='https://github.com/sheriftawfikabbas/pyphotonics',
    keywords=['DFT', 'Material science', 'Photoluminescence',
              'VASP'],
    install_requires=['scipy', 'numpy', 'pandas', 'matplotlib', 'pymatgen'],

)
