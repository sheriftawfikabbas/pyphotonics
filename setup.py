from distutils.core import setup
setup(
    name='pyphotonics',         # How you named your package folder (MyLib)
    packages=['pyphotonics'],   # Chose the same as "name"
    version='0.1.1',      # Start with a small number and increase it with every change you make
    # Chose a license from here: https://help.github.com/articles/licensing-a-repository
    license='gpl-3.0',
    description='The PyPhotonics python code is a post-processing code written entirely in python which takes as input the output files of the VASP and phonopy codes for a defect system, and calculates the Huang-Rhys factor and the PL lineshapes for that system.',   # Give a short description about your library
    author='Sherif Abdulkader Tawfik',                   # Type in your name
    author_email='sherif.tawfic@gmail.com',      # Type in your E-Mail
    # Provide either the link to your github or to your website
    url='https://github.com/sheriftawfikabbas/pyphotonics',
    # I explain this later on
    download_url='https://github.com/user/reponame/archive/v_01.tar.gz',
    keywords=['DFT', 'Material science', 'Photoluminescence',
              'VASP'],   # Keywords that define your package best
    install_requires=['scipy', 'numpy', 'pandas', 'matplotlib'],

)
