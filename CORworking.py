import main,plotters
import matplotlib.pyplot as plt
import numpy as np

resistances=main.linearResistance(10,0,10,0,1)
#plotters.PlotResistances(resistances)

growthRateArr = main.costOfSwitches(10,0.1,0.15)
#plotters.PlotResistances(growthRateArr)
model=main.GenModelCustGrowthRate(0.1,10, resistances,growthRateArr, 1, 400, [True]*7+[False]*7)
res = model.Run()
print(model.RunEig(True)[1])
model=main.GenModelCustGrowthRate(0,10, resistances,growthRateArr, 1, 400, [True])
resTwo = model.Run()

plotters.PlotTotalPop(res,"Total Population Over Time (Pulse Therapy)")
plotters.PlotTotalPop(resTwo,"Total Population Over Time (Continuous Therapy)")
#plotters.epiLinePlot(res,"Adaptive Therapy Subpopulation Size Over Time")
#plotters.epiLinePlot(resTwo,"Continuous Therapy Subpopulation Size Over Time")

print(np.sum(res)-np.sum(resTwo))

totArr = np.zeros(len(res))
for index in range(0,len(res)):
    totArr[index] = np.sum(res[index]) - np.sum(resTwo[index])
plt.plot(totArr)
plt.title("Total Population Difference (Adaptive minus Continuous)")
plt.xlabel("Days"); plt.ylabel("Population Difference")
plt.show()

plotters.EpiStackPlot(res)
plt.title("Phenotype Stack Plot (Adaptive Therapy)")
plt.xlabel("Days")
plt.ylabel("Population Size")
plt.show()

plotters.EpiStackPlot(resTwo)
plt.title("Phenotype Stack Plot (Continuous Therapy)")
plt.xlabel("Days")
plt.ylabel("Population Size")
plt.show()
# plotters.EpiStackPlot(res-resTwo)
# plt.title("Phenotype Stack Plot (Adaptive Therapy minus Continuous Therapy")
# plt.xlabel("Days")
# plt.ylabel("Population Size")
# plt.show()