import main
import numpy as np
import matplotlib.pyplot as plt

def objFun(treatmentMaxLen,holidayMaxLen,days,size):
    xs = np.arange(1,treatmentMaxLen+1,1)
    ys = np.arange(1,holidayMaxLen+1,1)
    resArr = np.zeros((len(xs),len(ys)))
    resistances = main.linearResistance(10,0,10,0,1)
    growthRateArr = main.costOfSwitches(10,0.1,0.15)
    for x,numOfZeros in enumerate(xs):
        for y,numOfOne in enumerate(ys):
            seq = [1]*numOfOne + [0]*numOfZeros; seqLen = len(seq)
            drugSeq = np.tile(seq,int(days/seqLen)+1)
            model=main.GenModelCustGrowthRate(.1,size, resistances,growthRateArr, 1, days, drugSeq)
            runSol = model.Run()
            resArr[x,y] = np.sum(runSol)
            if(x==0 or x==len(xs)-1) and (y==0 or y==len(ys)-1):
                print(str(x)+","+str(y)+" score:"+str(resArr[x,y]))

    return resArr.T,xs,ys

resThree,xs,ys = objFun(10,10,150,10)

resThree/= 150

plt.imshow(resThree)
minVal = [1000]
for i in range(0,len(resThree)):
    for j in range(0,len(resThree[i])):
        if resThree[i][j] < minVal[0]:
            minVal = [resThree[i][j],j+1,i+1]
print(minVal)
xAxis = np.arange(0, len(xs))
xs = [str(round(elem, 2)) for elem in xs]
plt.xticks(xAxis, xs)
plt.yticks(np.arange(0, len(ys)), ys)
plt.xlabel("Dosage On Time (days)")
plt.ylabel("Dosage Off Time (days)")
plt.title("Pulse Therapy Heat Map Measuring Tumor Size")
x = plt.colorbar()
x.ax.set_ylabel("Tumor Size")
plt.show()
