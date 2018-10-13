from scipy.integrate import odeint
import numpy as np
import main
import matplotlib.pyplot as plt

info = []
def nextSimulation(initialPop, mRate, k, tolerance, size, days, sequence):
    popArr = main.calcPopDistribution(size, initialPop)
    resistanceArr = main.assignResistance(size, [0, 1])
    growthRates = np.full(size,1)
    mutationRatesUp, mutationRatesDown = main.detMutationRate(mRate,size)
    results = odeint(main.calculateChange, popArr, np.arange(0, 100, 1), (growthRates, mutationRatesUp, mutationRatesDown, k, resistanceArr, main.periodicDrugDay(sequence, days), np.zeros(size)))

    totalPops = np.sum(results, axis=0); sumPop = np.sum(totalPops)
    mean, standardDev = main.getMeanAndStandardDeviation(popArr)
    info.append([sumPop, mean, standardDev])

    # For next iteration (essentially a recursion)
    tolerance = tolerance - 1
    if tolerance > 0:
        nextSimulation(initialPop, mRate, k, tolerance, size+1, days, sequence)
    return info


nSwitches = 100
res=nextSimulation(.01,.1,1,nSwitches,5,100,[False])
allTotals = np.zeros(nSwitches); allMeans = np.zeros(nSwitches); allStds = np.zeros(nSwitches);
for i in range(len(res)):
    allTotals[i] = res[i][0]
    allMeans[i] = res[i][1]
    allStds[i] = res[i][2]
# generate plots
plt.plot(allTotals); plt.title("Population Summation Graph"); plt.ylabel("Population Size"); plt.xlabel("Number of Switches"); plt.yticks(np.arange(0,150,30)); plt.show();
plt.plot(allMeans); plt.title("Weighted Mean Graph"); plt.ylabel("Mean"); plt.xlabel("Number of Switches"); plt.yticks(np.arange(0,1,.25)); plt.show()
plt.plot(allStds); plt.title("Standard Deviation Graph"); plt.ylabel("Deviation"); plt.xlabel("Number of Switches"); plt.show()