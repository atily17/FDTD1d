import matplotlib.pyplot as plt
import matplotlib.animation as anm
import numpy as np
import os

EYmax =  1
EYmin = -1
HYmax =  0.01
HYmin = -0.01

class SpaceGraph:
    def __init__(self, solver,oDirectory="graph"):
        self.fig = plt.figure(figsize = (5, 4))
        self.ax1 = self.fig.add_subplot(111)
        self.ax2 = self.ax1.twinx()
        self.iDirectory = solver.directory
        self.oDirectory = oDirectory
        self.padding = solver.padding
        self.N = len(solver.etimes)
        self.etimes = solver.etimes
        self.htimes = solver.htimes
        self.grid = solver.grid
        self.makeDirectory()

    def makeDirectory(self):
        os.makedirs(self.oDirectory, exist_ok=True)


    def plot(self,n,bShow=True,bSave=False, bCla=False, bClose=True):
        if bCla:
            self.ax1.cla()
            self.ax2.cla()
        filename = str(n).rjust(self.padding, "0")
        E = np.load(self.iDirectory + "\\E\\"+filename+".npy")
        H = np.load(self.iDirectory + "\\H\\"+filename+".npy")
        self.ax1.set_xlim([0, self.grid.domain])
        self.ax1.set_ylim([EYmin, EYmax])
        self.ax2.set_ylim([HYmin, HYmax])
        self.ax1.plot(self.grid.e_element,E, label="E",color="g")
        self.ax2.plot(self.grid.h_element,H, label="H",color="magenta")
        handler1, label1 = self.ax1.get_legend_handles_labels()
        handler2, label2 = self.ax2.get_legend_handles_labels()
        self.ax1.legend(handler1 + handler2, label1 + label2,  borderaxespad=0.)
        if bShow:
            self.fig.show()
        if bSave:
            self.fig.savefig(self.oDirectory+"\\EH"+filename+".png")
        if bClose:
            plt.close(self.fig)

    def update(self,i):
        self.plot(i, bShow=False,bCla=True, bClose=False)
        

    def animate(self):
        ani = anm.FuncAnimation(self.fig, self.update,  interval = 50, frames = self.N)
        ani.save(self.oDirectory+"\\animation.gif", writer = 'imagemagick')
        plt.close(self.fig)
        
    def load(self,i,eh):
        if eh == "E":
            f = np.load("E")
        elif eh == "H":
            f = np.load("H")
        return f


class FourierGraph:
    def __init__(self, solver,oDirectory="graph"):
        self.fig = plt.figure(figsize = (5, 4))
        self.ax1 = self.fig.add_subplot(111)
        self.iDirectory = solver.directory
        self.oDirectory = oDirectory
        self.padding = solver.padding
        self.grid = solver.grid
        self.freqs = solver.freqs
        self.fourierE = solver.fourierE
        self.fourierH = solver.fourierH
        self.makeDirectory()

    def makeDirectory(self):
        os.makedirs(self.oDirectory, exist_ok=True)

    def plotE(self, bShow=True, bSave=False):
        self.ax1.set_xlim((self.freqs[0], self.freqs[-1]))
        self.ax1.set_xscale("log")

        for key in self.fourierE:
            self.ax1.plot(self.freqs, np.abs(self.fourierE[key]), label=key)
        if bShow:
            plt.show()
        if bSave:
            self.fig.savefig(self.oDirectory+"\\fourierE.png")

    def plotH(self, bShow=True, bSave=False, bClose =True):
        self.ax1.set_xlim((self.freqs[0], self.freqs[-1]))
        self.ax1.set_xscale("log")

        for key in self.fourierH:
            self.ax1.plot(self.freqs, np.abs(self.fourierH[key]), label=key)
        if bShow:
            plt.show()
        if bSave:
            self.fig.savefig(self.oDirectory+"\\fourierH.png")
        if bClose:
            plt.close(self.fig)
