from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import main


def runSim(self, days):
	therapyArr = main.therapyListGenerator(days)
	singleCombination = [];
	popArr = []
	for i in range(0, len(therapyArr)):
		for elm in range(0, len(therapyArr[i]) - 1):
			singleCombination.append(therapyArr[i][elm])
		popArr.append(simulateDiffTreatments(self.popArr.size, self.popArr, self.k, self.growthRates, self.mRatesUp,
											 self.mRatesDown, self.resistances, singleCombination, days))
		therapyArr[i][len(therapyArr[i]) - 1] = popArr[i]
		singleCombination.clear()
	sortedTherapyArr = sorted(therapyArr, key=lambda x: x[len(therapyArr[0]) - 1])
	sortedTherapyArr = np.array(sortedTherapyArr);
	popArr = sorted(popArr)
	return sortedTherapyArr, popArr

#Purpose: calls the simulationDiffTreatment function for different treatment scedules
#Variable: sequence is the combination of drug therapy administration
def simulateDiffTreatments(size,popArr,k,growthRates,mutationRatesUp,mutationRatesDown,resistanceArr,sequence,days):
	resultArr = odeint(main.calculateChange,popArr,np.arange(0,days,1),(growthRates,mutationRatesUp,mutationRatesDown,k,resistanceArr,main.periodicDrugDay(sequence,days),np.zeros(size)))
	resultArr = np.sum(resultArr[-1, :])
	return resultArr

model=main.GenDefaultModel(0.05,15,1,100,[True])
results,vals=runSim(model,10)
plt.plot(vals)
plt.title("Simulating Different Treatment Schedules"); plt.xlabel("Treatment Rank"); plt.ylabel("Final Population Size")
plt.show()