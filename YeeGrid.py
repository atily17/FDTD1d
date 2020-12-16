import numpy as np
from scipy import constants

class Medium:
    def __init__(self):
        self.setMediumInit()
        self.sgm = 0
        self.left = 0
        self.right = 0

    def setRegion(self, left, right):
        self.left=left
        self.right=right

    def setMediumInit(self):
        self.eps = constants.epsilon0
        self.mu = constants.mu0

    def setMedium(self, reps, rmu, sgm):
        self.setMediumInit()
        self.eps = reps*self.eps
        self.mu = rmu*self.mu
        self.sgm = sgm

class YeeGrid:
    def __init__(self, directory="output"):
        self.domain = 0
        self.n = 0
        self.dz = 0
        self.dt = 0
        self.mediums = []
        self.e_element = np.array([])
        self.h_element = np.array([])

    def setGrid(self, d, N, ratio=0.99):
        self.domain = d
        self.n = N
        self.dz = self.domain/self.n
        self.e_element=np.linspace(0,self.domain, N+1)
        self.h_element=np.linspace(0+self.dz/2,self.domain-self.dz/2, self.n)
        c = constants.c
        self.dt = self.dz/c*ratio

    def addMedium(self, left, right, reps, rmu, sgm):
        m = Medium()
        m.setRegion(left, right)
        m.setMedium(reps, rmu, sgm)
        self.mediums.append(m)

    def getMediumsProperty(self):
        eps = np.full(len(self.e_element), constants.epsilon_0)
        mu = np.full(len(self.h_element), constants.mu_0)
        sgm = np.zeros(len(self.e_element))
        for medium in self.mediums:
            e_range = (medium.left < self.e_element) & (self.e_element < medium.right)
            h_range = (medium.left < self.h_element) & (self.h_element < medium.right)
            eps[e_range] = medium.eps
            sgm[e_range] = medium.sgm
            mu[h_range] = medium.mu

        properties = {"eps":eps, "sgm":sgm, "mu":mu}
        return properties

    def print(self):
        print("e:", self.e_element)
        print("h:", self.h_element)


