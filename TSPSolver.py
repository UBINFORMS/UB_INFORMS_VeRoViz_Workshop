from gurobipy import *

class CNode:
	def __init__(self, nodeID, nodeType, lat, lon):
		self.nodeID = nodeID
		self.nodeType = nodeType
		self.lat = lat
		self.lon = lon

def TSPSolver(nodes, tauDist, tauTime):
	n = len(nodes)

	TSP = Model('TSP')

	x = {}
	for i in range(n):
		for j in range(n):
			if i != j:
				x[i,j] = TSP.addVar(vtype=GRB.BINARY, obj = tauDist[i + 1, j + 1], name = 'x_{}_{}'.format(i,j))
				
	TSP.modelSense = GRB.MINIMIZE
	TSP.Params.lazyConstraints = 1
	TSP.update()

	for i in range(n):
		TSP.addConstr(quicksum(x[i,j] for j in range(n) if i != j) == 1, name='leave_{}'.format(i))
		TSP.addConstr(quicksum(x[j,i] for j in range(n) if i != j) == 1, name='enter_{}'.format(i))

	TSP._x = x
	TSP._n = n

	def subtourelim(model, where):
		if where == GRB.Callback.MIPSOL:
			x_sol = model.cbGetSolution(model._x)
			arcs = tuplelist((i,j) for i,j in model._x.keys() if x_sol[i,j] > 0.9)
			adjList = [[] for i in range(model._n)]
			for i, j in arcs:
				adjList[i].append(j)	
			found = [0 for i in range(model._n)]
			components = []
			for i in range(model._n):
				component = []
				queue = []
				if found[i] == 0:
					found[i] = 1
					component.append(i)
					queue.append(i)
					while queue:
						v = queue.pop(0)
						for u in adjList[v]:
							if found[u] == 0:
								found[u] = 1
								component.append(u)
								queue.append(u)
					components.append(component)
			for component in components:
				if len(component) < model._n:
					model.cbLazy(quicksum(x[i,j] for i in component for j in component if i != j) <= len(component) - 1)

	TSP.optimize(subtourelim)

	seq = []
	arcs = []
	if TSP.status == GRB.status.OPTIMAL:
		for i, j in x:
			if (x[i, j].x > 0.5):
				arcs.append([i, j])
		currentNode = 0
		currentTime = 0
		seq.append(currentNode + 1)
		while (len(arcs) > 0):
			for i in range(len(arcs)):
				if (arcs[i][0] == currentNode):
					currentNode = arcs[i][1]
					seq.append(currentNode + 1)
					arcs.pop(i)
					break

	return seq
