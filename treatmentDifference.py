import numpy as np
import matplotlib.pyplot as plt
import main,plotters

model=main.GenDefaultModel(0.2,10,1,100,[True]+[False]*10)
res = model.Run()
for i in range(0,len(res)):
    plotters.PlotHist(res[i])
    plt.title("Pulsed Therapy Population Distribution Over Time")
    plt.xlabel("Number of Switches On")
    plt.ylabel("Population Size")
    plt.ylim([0,.25])
    plt.show()
# iterations = len(res)
# trackArr = np.zeros([len(res),len(res[0])])
# for timeStep in range(0,iterations):
#     for subPop in range(0,len(res[timeStep])):
#         if timeStep > 0:
#             trackArr[timeStep][subPop] = res[timeStep][subPop]-res[timeStep-1][subPop]
#         else:
#             trackArr[timeStep][subPop] = 0
# changeArr = []
# for i in range(0,10):
#     for j in range(0,50):
#         changeArr.append(trackArr[j][i])
#     plt.plot(changeArr)
#     plt.title("Change of Population For Adaptive Therapy (switch " + str(i) + ")")
#     plt.xlabel("Time in Days")
#     plt.ylabel("Change in population")
#     plt.show()
#     del changeArr[:]