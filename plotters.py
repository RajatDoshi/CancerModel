import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.lines import Line2D

#Purpose: plot results with Rafael's special graphing code
#Variable: results is the vector that will be graphed
def EpiStackPlot(results):
	cmap=matplotlib.cm.get_cmap('inferno')
	mutColors=[cmap(x) for x in np.linspace(0,1,results.shape[1])]
	plt.stackplot(np.arange(0,results.shape[0]),results.T,zorder=1,colors=mutColors)

#Graph Histogram
def PlotHist(pops):
	plt.bar(np.arange(len(pops)),pops)

def PlotEigVals(results,xs,ys):
    print(xs)
    print(ys)
    plt.imshow(results)
    plt.xticks(np.arange(0,len(xs)),xs)
    plt.yticks(np.arange(0,len(ys)),ys)
    plt.colorbar()

def PlotResistances(resistances):
	plt.plot(resistances)
	plt.xlabel("# switches")
	plt.ylabel("resistance")
	plt.title("resistance plot")
	plt.show()

def PlotTotalPop(results,title):
	plt.plot(np.sum(results,axis=1))
	plt.title(title)
	plt.xlabel("Days")
	plt.ylabel("Total Population")
	plt.show()

def epiLinePlot(res,title):
	cmap=matplotlib.cm.get_cmap('coolwarm')
	mutColors=[cmap(x) for x in np.linspace(0,1,res.shape[1])]
	for i in range(res.shape[1]):
		plt.plot(res[:, i],color=mutColors[i])
	plt.title(title)
	plt.xlabel("Days")
	plt.ylabel("Subpopulation Size")
	plt.show()

def GenEpiLegend(nPops):
	cmap = matplotlib.cm.get_cmap('coolwarm')
	mutColors = [cmap(x) for x in np.linspace(0, 1, nPops)]
	labels = [str(x) + " switches on" for x in range(nPops)]
	GenLegend(mutColors,labels)

def GenLegend(colors,labels):

	handles=[Line2D([0], [0], color=colors[i], lw=4, label=labels[i]) for i in range(len(colors))]
	fig, ax = plt.subplots()
	ax.legend(handles=handles, loc='center')
	plt.show()



