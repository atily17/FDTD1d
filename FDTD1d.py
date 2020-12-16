import YeeGrid
import Solver
import Graph
import matplotlib.pyplot as plt
import numpy as np

N=1000

grid = YeeGrid.YeeGrid()
grid.setGrid(1, N)

solver = Solver.Solver(grid)
solver.setInit(1,100.0e-12)
solver.setBorder("absorb", "absorb")
solver.calc()

xe = grid.e_element
xh = grid.h_element

graph = Graph.SpaceGraph(solver)
graph.plot(0, bSave=True)
#graph.animate()

#graph = Graph.FourierGraph(solver)
#graph.plotE(bSave=True)
