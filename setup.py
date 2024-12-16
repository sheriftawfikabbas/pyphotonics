from setuptools import setup
from pathlib import Path

here = Path(__file__).resolve().parent
README = (here / "README.md").read_text(encoding="utf-8")
VERSION = (here / "pyphotonics" / "VERSION").read_text(encoding="utf-8").strip()

setup(
    name="pyphotonics",
    packages=["pyphotonics"],
    entry_points={
        "console_scripts": [
            "pyphotonics=pyphotonics.cli:execute_cli",
            "pyphotonics-incar=pyphotonics.cli:execute_incars",
        ],
    },
    include_package_data=True,
    version=VERSION,
    license="gpl-3.0",
    description="The PyPhotonics python code is a post-processing code written entirely in python which takes as input the output files of the VASP and phonopy codes for a defect system, and calculates the Huang-Rhys factor and the PL lineshapes for that system.",  # Give a short description about your library
    author="Sherif Abdulkader Tawfik",
    author_email="sherif.tawfic@gmail.com",
    long_description=README,
    long_description_content_type='text/markdown',
    url="https://github.com/sheriftawfikabbas/pyphotonics",
    keywords=["DFT", "Material science", "Photoluminescence", "VASP"],
    install_requires=[
        "scipy",
        "numpy",
        "pandas",
        "matplotlib",
        "pymatgen",
        "oganesson",
    ],
)
