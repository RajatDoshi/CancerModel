import main,plotters
import matplotlib.pyplot as plt
import numpy as np

resistances=main.linearResistance(6,0,4,0,1)
growthRateArr = main.costOfSwitches(6,0.5,0.5)
model=main.GenModelCustGrowthRate(0.01,6, resistances,growthRateArr, 1, 300, [True]*7+[False]*7)
res = model.Run()
print(model.RunEig(True)[1])
model=main.GenModelCustGrowthRate(0,6, resistances,growthRateArr, 1, 300, [True])
resTwo = model.Run()
print(np.sum(res)-np.sum(resTwo))

# totArr = np.zeros(len(res))
# for index in range(0,len(res)):
#     totArr[index] = np.sum(res[index]) - np.sum(resTwo[index])
# plt.plot(totArr)
# plt.title("Total Population Difference (Adaptive minus Continuous)")
# plt.xlabel("Days"); plt.ylabel("Population Difference")
# plt.show()

plt.plot(res)
plt.title("Phenotype Stack Plot (Adaptive Therapy)")
plt.xlabel("Days")
plt.ylabel("Population Size")
plt.show()

plt.plot(resTwo)
plt.title("Phenotype Stack Plot (Continuous Therapy)")
plt.xlabel("Days")
plt.ylabel("Population Size")
plt.show()
# plotters.EpiStackPlot(res-resTwo)
# plt.title("Phenotype Stack Plot (Adaptive Therapy minus Continuous Therapy")
# plt.xlabel("Days")
# plt.ylabel("Population Size")
# plt.show()