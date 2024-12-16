import argparse
import os
import random
import sys
from pathlib import Path
from typing import Optional
from .version import VERSION
from oganesson.io.vasp import Outcar


class CLI_Photoluminescence:
    def __init__(self, argv: Optional[str] = None) -> None:
        self.argv = argv or sys.argv[:]
        self.prog_name = Path(self.argv[0]).name
        self.formatter_class = argparse.RawDescriptionHelpFormatter
        self.epilog = ""
        self.version = VERSION

        default_locale = os.environ.get("LANG", "en_US").split(".")[0]
        # if default_locale not in AVAILABLE_LOCALES:
        #     default_locale = DEFAULT_LOCALE

        parser = argparse.ArgumentParser(
            prog=self.prog_name,
            description=f"{self.prog_name} version {self.version}",
            epilog=self.epilog,
            formatter_class=self.formatter_class,
        )

        parser.add_argument(
            "--version", action="version", version=f"%(prog)s {self.version}"
        )

        parser.add_argument("-v", "--verbose", action="store_true", help="verbosity")

        parser.add_argument(
            "--contcar_ground_state",
            "-cgs",
            type=str,
            help="directory of the ground state CONTCAR file",
        )

        parser.add_argument(
            "--phonopy_path",
            type=str,
            default="./phonopy/",
            help="directory to the phonopy band.yaml file",
        )

        parser.add_argument(
            "--num_modes",
            "-m",
            type=int,
            help="number of vibrational modes",
        )

        parser.add_argument(
            "--resolution",
            "-r",
            type=int,
            help="resolution of the PL line-shape",
        )

        parser.add_argument(
            "--method",
            "-M",
            type=str,
            help="method used to generate the vibrational modes",
        )

        arguments = parser.parse_args(self.argv[1:])

        from .photoluminescence import Photoluminescence

        photoluminescence = Photoluminescence(
            ground_state=arguments.ground_state,
            exceited_state=arguments.exceited_state,
            numModes=arguments.numModes,
            method=arguments.method,
            resolution=arguments.resolution,
            phonopy_path=arguments.phonopy_path,
        )

        photoluminescence.run()


class CLI_INCARs:
    def __init__(self, argv: Optional[str] = None) -> None:
        self.argv = argv or sys.argv[:]
        self.prog_name = Path(self.argv[0]).name
        self.formatter_class = argparse.RawDescriptionHelpFormatter
        self.epilog = ""
        self.version = VERSION

        default_locale = os.environ.get("LANG", "en_US").split(".")[0]
        # if default_locale not in AVAILABLE_LOCALES:
        #     default_locale = DEFAULT_LOCALE

        parser = argparse.ArgumentParser(
            prog=self.prog_name,
            description=f"{self.prog_name} version {self.version}",
            epilog=self.epilog,
            formatter_class=self.formatter_class,
        )

        parser.add_argument(
            "--version", action="version", version=f"%(prog)s {self.version}"
        )

        parser.add_argument("-v", "--verbose", action="store_true", help="verbosity")

        parser.add_argument(
            "-f",
            "--vasp_folder",
            type=str,
            default="./",
            help="the location of the VASP folder",
        )

        arguments = parser.parse_args(self.argv[1:])

        outcar = Outcar(arguments.vasp_folder, "OUTCAR")
        u, d = outcar.get_ferwe_ferdo()
        print("FERWE =", u)
        print("FERDO =", d)


def execute_cli(argv: Optional[str] = None) -> None:
    cli = CLI_Photoluminescence()


def execute_incars(argv: Optional[str] = None) -> None:
    cli = CLI_INCARs()


if __name__ == "__main__":
    execute_cli()
