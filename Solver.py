from enum import Enum
import numpy as np
from scipy import constants
import os

class BoundaryCond(Enum):
    ABSORB = 1
    DIRICHLET = 2

class Pulse:
    def __init__(self):
        self.w = 100e-9
        self.threshold = 4
        self.E = 1

    def calc(self,grid):
        c = constants.c
        mu0 = constants.mu_0
        eps0 = constants.epsilon_0
        impedance = np.sqrt(mu0/eps0)
        z0 = c*self.w*self.threshold+2*grid.dz
        taues = -(grid.e_element-z0)/c
        tauhs = -grid.dt/2 - (grid.h_element-z0)/c

        return {"E": self.E*self._gaussian(taues, self.w), "H":self.E/impedance*self._gaussian(tauhs, self.w)}

    def _gaussian(self, taus, w):
        threshold = self.threshold*w
        flter = (-threshold < taus) & (taus < threshold)
        p = np.zeros(len(taus))
        p[flter] = np.exp(-taus[flter]**2/(2*w**2))
        return p


class Solver:
    def __init__(self, grid, directory="output"):
        self.borderL = BoundaryCond.ABSORB
        self.borderR = BoundaryCond.ABSORB
        self.grid = grid
        self.E = np.zeros(len(self.grid.e_element))
        self.H = np.zeros(len(self.grid.h_element))
        self.NTime = 10000
        self.properties = self.grid.getMediumsProperty()
        self.time = 0
        self.etimes = [0]
        self.htimes = [-self.grid.dt/2]
        self.observes = []
        self.fourierE = {}
        self.fourierH = {}
        self.pulse=Pulse()
        self.directory = directory
        self.padding = 6
        self.makeDirectory()
        self.setFourier(10e6, 10e9, 90)
        self.addObserve(0.5)
        self.check = 10

    def makeDirectory(self):
        os.makedirs(self.directory, exist_ok=True)
        os.makedirs(self.directory + "\\E", exist_ok=True)
        os.makedirs(self.directory + "\\H", exist_ok=True)

    def save(self, n):
        filename = str(n).rjust(self.padding, "0")
        np.save(self.directory + "\\E\\"+filename, self.E)
        np.save(self.directory + "\\H\\"+filename, self.H)

    def setInit(self, E, w):
        self.pulse.E = E
        self.pulse.w = w

    def setBorder(self, l, r):
        if l == "absorb":
            self.borderL = BoundaryCond.ABSORB
        elif l == "dirichlet":
            self.borderL = BoundaryCond.DIRICHLET
        if r == "absorb":
            self.borderR = BoundaryCond.ABSORB
        elif r == "dirichlet":
            self.borderR = BoundaryCond.DIRICHLET

    def setFourier(self, fmin, fmax, n):
        fmin = np.log10(fmin)
        fmax = np.log10(fmax)
        self.freqs=np.logspace(fmin, fmax, num=n+1)
        
    def addObserve(self, z, name=""):
        eidx = (np.abs(self.grid.e_element - z)).argmin()
        hidx = (np.abs(self.grid.h_element - z)).argmin()
        if name=="":
            name="observe"+str(len(self.observes))
        idx = {"E": eidx, "H": hidx, "name": name}
        self.observes.append(idx)
        self.fourierE[name] = np.zeros(len(self.freqs), dtype=np.complex128)
        self.fourierH[name] = np.zeros(len(self.freqs), dtype=np.complex128)

    def calc(self):
        self.calcInit()
        self.save(0)
        self.time = 0
        self.calcFourier()
        for t in range(1, self.NTime+1):
            self.next(t)
            self.save(t)            
            if t % 10 == 0:
                print("iter="+str(t)+", ", end="")
                if self.judgeConvergence(bPrint=True):
                    break

    def next(self, t):
        preE = self.E.copy()
        preH = self.H.copy()
        self.calcH()
        self.time += self.grid.dt/2
        self.htimes.append(self.time)
        self.calcFourier("H")
        self.calcE()
        self.time += self.grid.dt/2
        self.etimes.append(self.time)
        self.calcBorder(preE)
        self.calcFourier("E")
            
    def calcInit(self):
        p = self.pulse.calc(self.grid)
        self.E = p["E"]
        self.H = p["H"]

    def calcH(self):
        mu = self.properties["mu"]
        dt = self.grid.dt
        dz = self.grid.dz
        c1 = 1/mu*dt/dz
        self.H = self.H - c1 * (self.E[1:] - self.E[:len(self.E)-1])

    def calcE(self):
        I=len(self.E)
        eps = self.properties["eps"][1:I-1]
        sgm = self.properties["sgm"][1:I-1]
        dt = self.grid.dt
        dz = self.grid.dz
        c1 = (2*eps-sgm*dt)/(2*eps+sgm*dt)
        c2 = 2*dt/((2*eps+sgm*dt)*dz)
        self.E[1:I-1] = c1*self.E[1:I-1] - c2*(self.H[1:]-self.H[:len(self.H)-1])

    def calcBorder(self, preE):
        c = constants.c
        dz = self.grid.dz
        dt = self.grid.dt
        if self.borderL == BoundaryCond.ABSORB:
            c1 = (c*dt - dz)/(c*dt+dz)
            self.E[0] = preE[1] +  c1 * (self.E[1] - preE[0])
        elif self.borderL == BoundaryCond.DIRICHLET:
            self.E[0] = 0

        I = len(self.E) - 1
        if self.borderR == BoundaryCond.ABSORB:
            c1 = (c*dt - dz)/(c*dt+dz)
            self.E[I] = preE[I-1] +  c1 * (self.E[I-1] - preE[I])
        elif self.borderR == BoundaryCond.DIRICHLET:
            self.E[I] = 0

    def calcFourier(self, eh=""):
        j = 1j
        f = np.exp(-j*self.time*2*np.pi*self.freqs)
        for observe in self.observes:
            if eh != "H":
                self.fourierE[observe["name"]] += self.E[observe["E"]]*f
            if eh != "E":
                self.fourierH[observe["name"]] += self.H[observe["H"]]*f

    def judgeConvergence(self, cc = 1e-3, bPrint=False):
        emax = self.pulse.E
        eave = np.linalg.norm(self.E, 1)/len(self.E)
        conv = eave/emax
        bConvergence = conv < cc
        if bPrint:
            print("time="+str(self.time)+", max/norm="+str(conv))
        return bConvergence

    def print(self):
        print("e:", self.E)
        print("h:", self.H)