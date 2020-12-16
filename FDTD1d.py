import YeeGrid
import Solver
import Graph

N=500

grid = YeeGrid.YeeGrid()
grid.setGrid(1, N)
grid.addMedium(0.4, 0.6, 3, 1, 0)

solver = Solver.Solver(grid)
solver.setInit(1,100.0e-12)
solver.setBorder("absorb", "absorb")
solver.calc()

graph = Graph.SpaceGraph(solver)
graph.animate(intervl=50)

graph = Graph.FourierGraph(solver)
graph.plotE(bSave=True)
