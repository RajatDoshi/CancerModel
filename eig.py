import main
import numpy as np
import matplotlib.pyplot as plt
import plotters

def eigVary(initSize,endSize,mRateStart,mRateEnd,step):
    sizeArr = np.arange(initSize, endSize + 1, 1)
    mRateArr = np.arange(mRateStart, mRateEnd + step, step)
    results = np.zeros((len(sizeArr), len(mRateArr)))
    for x, mRate in enumerate(mRateArr):
        for y, size in enumerate(sizeArr):
            resistances = main.expResAssign(size,1.0001,0,size-1)
            model=main.GenModelCustRes(mRate*size,size,resistances,1,100,[True])
            res = model.RunEig(True)
            results[y,x] = res[1]
    return results, mRateArr, sizeArr

#res,xs,ys = eigVary(2,15,.01,.2,.01)
# plt.imshow(res)
# xAxis = np.arange(0, len(xs))
# xs = [str(round(elem, 2)) for elem in xs]
# xs = [elem.strip("0") for elem in xs]
# plt.xticks(xAxis, xs)
# plt.yticks(np.arange(0, len(ys)), ys)
# plt.xlabel("Mutation Rate")
# plt.ylabel("Number of Switches")
# plt.title("Steady State Growth Rates")
# plt.colorbar()
# plt.show()

resistances=main.linearResistance(10,0,10,0,1)
model = main.GenModelCustRes(.2, 10, resistances, 1, 100, [False])
res = model.RunEig(False)
plotters.PlotHist(res[2])
plt.xlabel("Number of Switches On"); plt.xticks(np.arange(0,len(res[2]),1))
plt.ylabel("Relative Distribution"); plt.title("Population Distribution (No Drug Applied)")
plt.show()