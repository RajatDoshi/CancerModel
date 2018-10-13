import main,plotters
import matplotlib.pyplot as plt
import numpy as np

resistances=main.linearResistance(10,0,10,0,1)
# plotters.PlotResistances(resistances)

growthRateArr = main.costOfSwitches(10,0.5,0.5)
# plotters.PlotResistances(growthRateArr)
model=main.GenModelCustGrowthRate(0.1,10, resistances,growthRateArr, 1, 400, [True]*7+[False]*7)
res = model.Run()
print(model.RunEig(True)[1])
model=main.GenModelCustGrowthRate(0,10, resistances,growthRateArr, 1, 400, [True])
resTwo = model.Run()

print(np.sum(res)-np.sum(resTwo))

totArr = np.zeros(len(res))
for index in range(0,len(res)):
    totArr[index] = np.sum(res[index]) - np.sum(resTwo[index])
plt.plot(totArr)
plt.title("Tumor Size Difference (Pulsed Minus Continuous Therapy")
plt.xlabel("Days"); plt.ylabel("Population Difference")
plt.show()

plotters.PlotTotalPop(res,"Total Population Over Time (Pulsed Therapy)")
plotters.PlotTotalPop(resTwo,"Total Population Over Time (Continuous Therapy)")

plotters.epiLinePlot(res,"Pulsed Therapy Subpopulation Size Over Time")
plotters.epiLinePlot(resTwo,"Continuous Therapy Subpopulation Size Over Time")
plotters.GenEpiLegend(res.shape[1])

labels=[str(x)+" switches on" for x in range(res.shape[1])]
for i in range(res.shape[1]):
    plt.plot(res[:,i],label=labels[i])
plt.legend(labels=labels)
plt.title("Pulsed Therapy Subpopulation Change over Time")
plt.xlabel("Days")
plt.ylabel("Subpopulation Size")
plt.show()

plt.plot(resTwo)
plt.title("Continuous Therapy Subpopulation Change over Time")
plt.xlabel("Days")
plt.ylabel("Subpopulation Size")
plt.show()
plotters.EpiStackPlot(res-resTwo)
plt.title("Phenotype Stack Plot (Pulsed Therapy minus Continuous Therapy")
plt.xlabel("Days")
plt.ylabel("Population Size")
plt.show()