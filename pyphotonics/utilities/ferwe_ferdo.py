from oganesson.io.vasp import Outcar


class FerweFerdo:
    def __init__(self, path="./"):
        incarf = open(path + "INCAR")
        incar = incarf.read()
        incarf.close()
        outcar = Outcar(path)
