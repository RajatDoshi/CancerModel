from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import main

#Purpose: test varying parameters (resistance and mRate)
#Variable: isMin is a boolean for whether you are changing the min or the max resistance
#Input resistance values that align to switch number startng at 1
def varyParameters(initialMutRate,endMutRate,resistanceArr,size,growthRates,isDrug,isMin):
    step = .05
    xs = np.arange(initialMutRate,endMutRate+step,step)
    ys = np.arange(resistanceArr[0],resistanceArr[1]+1,1)
    results=np.zeros((len(ys),len(xs)))
    for x,mRate in enumerate(xs):
        for y,resVal in enumerate(ys):
            mutationRateUp, mutationRateDown = main.detMutationRate(mRate,size)
            if isMin:
                resistances = main.resAssign(size,[0,1],resVal,size-1)
            else:
                resistances = main.resAssign(size, [0, 1], 0, resVal)
                print(resistances)
            eigVal = main.getEig(growthRates, resistances, mutationRateUp, mutationRateDown, isDrug)
            results[y,x]=eigVal[1]
    return results,xs,ys

#Purpose: vary mRate and resistance and plot eig
def varyExpMRate(initialMutRate,endMutRate,resistanceArr,n,size,growthRates,isDrug,isMin):
    step = .01; xs = np.arange(initialMutRate,endMutRate,step); ys = np.arange(resistanceArr[0],resistanceArr[1],1)
    results=np.zeros((len(ys),len(xs)))
    for x,mRate in enumerate(xs):
        mutationRateUp, mutationRateDown = main.detMutationRate(mRate, size)
        for y,resVal in enumerate(ys):
            if isMin:
                resistances = main.expResAssign(size,n,resVal,size)
            else:
                resistances = main.expResAssign(size,n,0,resVal)
            eigVal = main.getEig(growthRates, resistances, mutationRateUp, mutationRateDown, isDrug)
            results[y,x]=eigVal[1]
    return results,xs,ys

#Purpose: vary mutation rate and nSwitches and total plot eig
def varyParametersSize(initialMutRate,endMutRate,maxSize,isDrug):
    xs = np.arange(initialMutRate,endMutRate+.01,.01)
    ys = np.arange(2,maxSize+1,1)
    results=np.zeros((len(ys),len(xs)))
    for x,mRate in enumerate(xs):
        for y,size in enumerate(ys):
            mutationRateUp, mutationRateDown = main.detMutationRate(mRate*size,size)
            growthRates = np.full(size,.5)
            resistanceVal = main.expResAssign(size,1.15,0,size-1)
            #plotters.PlotHist(resistanceVal)
            #plt.xlabel("State"); plt.ylabel("Resistance Values"); plt.title("Resistance Histogram")
            #plt.show()
            eigVal = main.getEig(growthRates, resistanceVal, mutationRateUp, mutationRateDown,isDrug)
            results[y,x]=eigVal[1]
    return results,xs,ys

#Purpose: vary mRate and nSwitches and plot totalPopSize
def varyParametersPop(initialMutRate,endMutRate,maxSize,k,drugWeeks,timesteps,initPopSize):
    xs = np.arange(initialMutRate,endMutRate,.01)
    ys = range(2,maxSize+1)
    resultsArr = np.zeros((len(ys),len(xs)))
    for x,mRate in enumerate(xs):
        for y,size in enumerate(ys):
            mutationRateUp, mutationRateDown = main.detMutationRate(mRate*size,size)
            popArr = main.calcPopDistribution(size,initPopSize)
            resistances = main.assignResistance(size, [0, 1])
            change = np.zeros(size); growthRates = np.full(size,1)
            results = odeint(main.calculateChange, popArr, timesteps, (growthRates, mutationRateUp, mutationRateDown, k, resistances, drugWeeks, change))
            resultsArr[y,x] = np.sum(results[:][-1])
    return resultsArr,xs,ys

def RunVary(startMutRate, endMutRate, maxSize, isDrug):
    return varyParametersSize(startMutRate, endMutRate, maxSize, isDrug)

def RunPop(self, startMutRate, endMutRate, maxSize, initPopSize):
    return varyParametersPop(startMutRate, endMutRate, maxSize, self.k, self.drugWeeks, self.timesteps, initPopSize)

def RunExpResMRate(self, startMutRate, endMutRate, resistanceBounds, n, isDrug, isMin):
    return varyExpMRate(startMutRate, endMutRate, resistanceBounds, n, self.popArr.size, self.growthRates, isDrug,
                            isMin)

def RunMin(self, startMutRate, endMutRate, resistanceMaxBounds, isDrug, isMin):
    return varyParameters(startMutRate, endMutRate, resistanceMaxBounds, self.popArr.size, self.growthRates, isDrug,
                              isMin)

model=main.GenDefaultModel(0.1,15,1,50,[False])
resThree,xs,ys = RunVary(.01,.2,15,True)
plt.imshow(resThree)
xAxis = np.arange(0, len(xs))
xs = [str(round(elem, 2)) for elem in xs]
xs = [elem.strip("0") for elem in xs]
plt.xticks(xAxis, xs)
plt.yticks(np.arange(0, len(ys)), ys)
plt.xlabel("Transition Rate")
plt.ylabel("Number of Switches")
plt.title("Steady State Growth Rate under Constant Drug Application")
plt.colorbar()
plt.show()