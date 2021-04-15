import numpy as np
from numpy import fft
from XYZ import XYZ
import sys
import matplotlib.pyplot as plt
import cmath


class Photoluminescence:

    def vasp_read_modes(self):
        return 0

    def phonopy_read_modes(self):
        modes = np.zeros((self.numModes, self.numAtoms, 3))

        try:
            band = open(self.path + "band.yaml", 'r')
        except OSError:
            print("Could not open/read file: band.yaml")
            sys.exit()

        for line in band:
            if "  band:" in line:
                break

        for i in range(self.numModes):
            band.readline()
            band.readline()
            band.readline()

            for a in range(self.numAtoms):
                line = band.readline()

                line = band.readline().replace(",", "")
                parts = line.strip().split()
                modes[i][a][0] = float(parts[2])

                line = band.readline().replace(",", "")
                parts = line.strip().split()
                modes[i][a][1] = float(parts[2])

                line = band.readline().replace(",", "")
                parts = line.strip().split()
                modes[i][a][2] = float(parts[2])

        band.close()

        return modes

    def phonopy_read_frequencies(self):
        frequencies = np.zeros(self.numModes)
        try:
            band = open(self.path + "band.yaml", 'r')
        except OSError:
            print("Could not open/read file: band.yaml")
            sys.exit()

        for line in band:
            if "  band:" in line:
                break

        for i in range(self.numModes):
            band.readline()
            line = band.readline()

            parts = line.strip().split()
            frequencies[i] = float(parts[1])

            line = band.readline()

            for a in range(self.numAtoms):
                band.readline()
                band.readline()
                band.readline()
                band.readline()

        band.close()
        return frequencies

    def vasp_read_frequencies(self):
        return 0

    def get_S_omega(self, omega, sigma):
        sum = 0
        for k in range(len(self.S)):
            sum += self.S[k] * self.gaussian(omega, self.frequencies[k], sigma)
        return sum

    def gaussian(self, omega, omega_k, sigma):
        return 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(omega - omega_k) * (omega - omega_k) / sigma / sigma / 2)

    def write_S(self, file_name):
        f = open(file_name, 'w')
        for i in range(len(self.S_omega)):
            # f.write(str(self.omega_set[i]) + "\t" + str(self.S_omega[i])+'\n')
            f.write(str(self.S_omega[i])+'\n')
        f.close()

    def PL(self, gamma, SHR, EZPL):
        Gt = []
        I = []

        r = 1/self.resolution
        St = fft.ifft(self.S_omega)
        St = fft.ifftshift(St)
        G = np.exp(2*np.pi*St-SHR)

        for i in range(len(G)):
            t = r*(i-len(G)/2)
            Gt += [G[i]*np.exp(-gamma*np.abs(t))]

        A = fft.fft(Gt)

        # Now, shift the ZPL peak to the EZPL energy value
        tA = A.copy()
        for i in range(len(A)):
            A[(int(EZPL*self.resolution)-i) % len(A)] = tA[i]

        for i in range(len(A)):
            I += [A[i]*((i)*r)**3]

        return A, np.array(I)

    def __init__(self, path, xyz_g, xyz_e, numModes, method, m, resolution):
        self.resolution = resolution
        self.numModes = numModes
        self.path = path
        self.m = m
        self.g = XYZ(xyz_g)
        self.e = XYZ(xyz_e)
        self.numAtoms = self.g.N
        R_g = np.zeros((self.numAtoms, 3))
        R_e = np.zeros((self.numAtoms, 3))
        self.method = method
        self.m = m

        for i in range(self.numAtoms):
            R_g[i][0] = self.g.x[i]
            R_g[i][1] = self.g.y[i]
            R_g[i][2] = self.g.z[i]
            R_e[i][0] = self.e.x[i]
            R_e[i][1] = self.e.y[i]
            R_e[i][2] = self.e.z[i]

        if "phonopy" in method:
            r = self.phonopy_read_modes()
            self.frequencies = self.phonopy_read_frequencies()
        else:
            r = self.vasp_read_modes()
            self.frequencies = self.vasp_read_frequencies()

        self.HuangRhyes = 0
        self.Delta_R = 0
        self.Delta_Q = 0
        self.IPR = []
        self.q = []
        self.S = []

        for i in range(numModes):
            q_i = 0
            IPR_i = 0
            participation = 0
            if method == "vasp":
                self.frequencies[i] = self.frequencies[i] / 1000
            elif method == "phonopy":
                self.frequencies[i] = self.frequencies[i] * \
                    0.004135665538536  # THz
            elif method == "phonopy-siesta":
                self.frequencies[i] = self.frequencies[i] * \
                    0.004135665538536 * 0.727445665  # THz

            if self.frequencies[i] < 0:
                self.frequencies[i] = 0

            max_Delta_r = 0
            for a in range(self.numAtoms):
                # Normalize r:
                participation = r[i][a][0] * r[i][a][0] + \
                    r[i][a][1] * r[i][a][1] + r[i][a][2] * r[i][a][2]
                IPR_i += participation**2

                for coord in range(3):
                    q_i += np.sqrt(m[a]) * (R_e[a][coord] -
                                            R_g[a][coord]) * r[i][a][coord] * 1e-10
                    if np.abs(r[i][a][coord]) > max_Delta_r:
                        max_Delta_r = np.abs(r[i][a][coord])

            IPR_i = 1.0 / IPR_i
            S_i = self.frequencies[i] * q_i**2 / 2 * 1.0 / \
                (1.0545718e-34 * 6.582119514e-16)

            self.IPR += [IPR_i]
            self.q += [q_i]
            self.S += [S_i]
            self.HuangRhyes += S_i

        for a in range(self.numAtoms):
            for coord in range(3):
                self.Delta_R += (R_e[a][coord] - R_g[a][coord])**2
                self.Delta_Q += (R_e[a][coord] - R_g[a][coord])**2 * m[a]

        self.Delta_R = self.Delta_R**0.5

        self.Delta_Q = (self.Delta_Q / 1.660539040e-27) ** 0.5

        self.max_energy = 5

        self.omega_set = np.linspace(
            0, self.max_energy, self.max_energy*self.resolution)
        self.S_omega = [self.get_S_omega(o, 6e-3) for o in self.omega_set]
        

    def print_table(self):
        for i in range(self.numModes):
            print("IPR\t", i, "\tSk\t", self.S[i], "\tenergy\t",
                  self.frequencies[i], "\t=\t", self.IPR[i], "\twith localization ratio beta =\t", 64 / self.IPR[i])
