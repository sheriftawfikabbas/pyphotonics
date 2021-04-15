import scipy.linalg as scl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

hbar = 1
m = 1
a = 20


class Schrodinger:
    def __init__(self, QE_data_csv):
        self.N = 2000
        self.QE_data = pd.read_csv(QE_data_csv, header=0, index_col=None)

        self.Q = self.QE_data.Q
        self.E = self.QE_data.E
        
        coeff = np.polynomial.polynomial.polyfit(
            self.QE_data.Q, self.QE_data.E, 6)

        Qmin = self.QE_data.Q.min()
        Qmax = self.QE_data.Q.max()

        roots=np.polynomial.polynomial.polyroots((coeff[1],2*coeff[2],3*coeff[3],4*coeff[4],5*coeff[5],6*coeff[6]))
        Qroots = []
        for r in roots:
            if np.imag(r) == 0 and (r > Qmin and r<Qmax):
                Qroots += r

        V_at_roots = [self.potential(coeff, xi) for xi in Qroots]

        self.Q0 = Qroots[V_at_roots.index(min(V_at_roots))]
        self.E0 = min(V_at_roots)
        self.QE_data.E = self.QE_data.E - self.E0

        x = np.linspace(-a/2., a/2., self.N)
        h = x[1]-x[0]
        V = [self.potential(coeff, xi) for xi in x]

        Mdd = 1./(h*h)*(np.diag(np.ones(self.N-1), -1) - 2 *
                        np.diag(np.ones(self.N), 0) + np.diag(np.ones(self.N-1), 1))
        H = -(hbar*hbar)/(2.0*m)*Mdd + np.diag(V)
        self.epsilon, psiT = np.linalg.eigh(H)
        self.chi = np.transpose(psiT)

        


        plt.figure(figsize=(10, 7))
        for i in range(5):
            # Flip the wavefunctions if it is negative at large x, so plots are more consistent.
            if self.chi[i][self.N-10] < 0:
                plt.plot(x, -self.chi[i]/np.sqrt(h),
                         label="$E_{}$={:>8.3f}".format(i, self.epsilon[i]))
            else:
                plt.plot(x, self.chi[i]/np.sqrt(h),
                         label="$E_{}$={:>8.3f}".format(i, self.epsilon[i]))
            plt.title("Solutions to the vibrational equation")
        plt.legend()
        plt.savefig("SE_WaveFunctions", bbox_inches='tight')


    def potential(self, coeff, x):
        return coeff[0]+x*coeff[1]+x**2*coeff[2]+x**3*coeff[3]+x**4*coeff[4]+x**5*coeff[5]+x**6*coeff[6]
