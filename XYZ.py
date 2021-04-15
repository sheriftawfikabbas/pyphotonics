class XYZ:
    def __init__(self, file_name):

        self.atom = []
        self.x = []
        self.y = []
        self.z = []
        

        f = open(file_name,'r')
        self.N = int(f.readline())
        f.readline()
        for i in range(self.N):
            line = f.readline().strip().split()
            self.atom += [line[0]]
            self.x += [float(line[1])]
            self.y += [float(line[2])]
            self.z += [float(line[3])]
