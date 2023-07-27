import numpy as np

class XYZ:
    def __init__(self, file_name):

        self.atom = []
        self.coordinates = []

        f = open(file_name,'r')
        self.N = int(f.readline())
        f.readline()
        for i in range(self.N):
            line = f.readline().strip().split()
            self.atom += [line[0]]
            
            c = [float(line[1]), float(line[2]), float(line[3])]
            self.coordinates += [c]
        
        self.coordinates = np.array(self.coordinates)
    def __len__(self):
        return len(self.atom)