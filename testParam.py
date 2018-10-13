from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import main

def tumorSizeCalcExp(sequence,mRate,size,reachMaxRes):
    days = len(sequence);k = 1; growthRates = np.full(size,.5)
    resistances = main.expResAssign(size,1.001,0,reachMaxRes)
    popArr = main.calcPopDistribution(int(size),.01); mutationRatesUp, mutationRatesDown = main.detMutationRate(mRate,size)
    resultArr = odeint(main.calculateChange, popArr, np.arange(0, days, 1), (growthRates, mutationRatesUp, mutationRatesDown, k, resistances, main.periodicDrugDay(sequence, days),np.zeros(popArr.size)))
    return np.sum(resultArr)

def tumorSizeCalcLinear(sequence,mRate,size):
    days = len(sequence);k = 1; growthRates = np.full(size,.5)
    resistances = main.assignResistance(size,[0,1])
    popArr = main.calcPopDistribution(int(size),.01); mutationRatesUp, mutationRatesDown = main.detMutationRate(mRate,size,True,False)
    resultArr = odeint(main.calculateChange, popArr, np.arange(0, days, 1), (growthRates, mutationRatesUp, mutationRatesDown, k, resistances, main.periodicDrugDay(sequence, days),np.zeros(popArr.size)))
    return np.sum(resultArr[50,:])

def tumorSizeCalcChangingMaxRes(sequence,mRate,size,maxResSwitchNumb,n):
    days = len(sequence); k = 1; growthRates = np.full(size,.5)
    resistances = main.expResAssign(size,n,0,maxResSwitchNumb)
    popArr = main.calcPopDistribution(int(size),.01); mutationRatesUp, mutationRatesDown = main.detMutationRate(mRate,size)
    resultArr = odeint(main.calculateChange, popArr, np.arange(0, days, 1), (growthRates, mutationRatesUp, mutationRatesDown, k, resistances, main.periodicDrugDay(sequence, days),np.zeros(popArr.size)))
    return np.sum(resultArr)

def costOfResFunction(sequence,mRate,size,maxResSwitchNumb,n):
    days = len(sequence); k = 1; resistances = main.expResAssign(size,n,0,maxResSwitchNumb)
    popArr = main.calcPopDistribution(int(size),1); mutationRatesUp, mutationRatesDown = main.detMutationRate(mRate,size)
    growthRate = main.costOfRes(.25,resistances)
    resultArr = odeint(main.calculateChange, popArr, np.arange(0, days, 1), (growthRate, mutationRatesUp, mutationRatesDown, k, resistances, main.periodicDrugDay(sequence, days),np.zeros(popArr.size)))
    return np.sum(resultArr)

#Purpose: Test different combinations of mRate and size parameters
#Variables: adaptiveTherapyBoolean is True if sequence is adaptive therapy but false if it is continous therapy
def varyCostRes(mRateStart,mRateEnd,step,initSize,endSize,days,adaptiveTherapyBoolean):
    sizeArr = np.arange(initSize,endSize+1,1)
    mRateArr = np.arange(mRateStart,mRateEnd+step,step)
    results=np.zeros((len(sizeArr),len(mRateArr)))
    for x,mRate in enumerate(mRateArr):
        for y,size in enumerate(sizeArr):
            if adaptiveTherapyBoolean:
                sequence = np.tile([True,False],(int(days/2)+1))
                results[y,x]=costOfResFunction(sequence,mRate,size,size-1,1.001)
            else:
                sequence = np.full(days,True)
                results[y,x]=costOfResFunction(sequence, mRate, size,size-1,1.001)
    return results,mRateArr,sizeArr

#Purpose: Test different combinations of mRate and size parameters
#Variables: adaptiveTherapyBoolean is True if sequence is adaptive therapy but false if it is continous therapy
def varyMRateSizeLinear(mRateStart,mRateEnd,step,initSize,endSize,days,adaptiveTherapyBoolean):
    sizeArr = np.arange(initSize,endSize+1,1)
    mRateArr = np.arange(mRateStart,mRateEnd+step,step)
    results=np.zeros((len(sizeArr),len(mRateArr)))
    for x,mRate in enumerate(mRateArr):
        for y,size in enumerate(sizeArr):
            if adaptiveTherapyBoolean:
                sequence = np.tile([True,False],int(days/2)+1)
                results[y,x]=tumorSizeCalcLinear(sequence,mRate,size)
            else:
                sequence = np.full(days,True)
                results[y,x]=tumorSizeCalcLinear(sequence, mRate, size)
    return results,mRateArr,sizeArr

#Purpose: Test different combinations of mRate and size parameters
#Variables: adaptiveTherapyBoolean is True if sequence is adaptive therapy but false if it is continous therapy
def varyMRateSizeExp(mRateStart,mRateEnd,step,initSize,endSize,days,adaptiveTherapyBoolean):
    sizeArr = np.arange(initSize,endSize+1,1)
    mRateArr = np.arange(mRateStart,mRateEnd+step,step)
    results=np.zeros((len(sizeArr),len(mRateArr)))
    for x,mRate in enumerate(mRateArr):
        for y,size in enumerate(sizeArr):
            if adaptiveTherapyBoolean:
                sequence = np.tile([True,False],int(days/2)+1)
                results[y,x]=tumorSizeCalcExp(sequence,mRate,size,size-1)
            else:
                sequence = np.full(days,True)
                results[y,x]=tumorSizeCalcExp(sequence, mRate, size,size-1)
    return results,mRateArr,sizeArr

#Purpose: Test different combinations of mRate and size parameters
#Variables: adaptiveTherapyBoolean is True if sequence is adaptive therapy but false if it is continous therapy
def varyMaxResStart(mRateStart,mRateEnd,step,n,initSize,endSize,days,numbOfSwitches,adaptiveTherapyBoolean):
    sizeArr = np.arange(initSize,endSize+1,1)
    mRateArr = np.arange(mRateStart,mRateEnd+step,step)
    results=np.zeros((len(sizeArr),len(mRateArr)))
    for x,mRate in enumerate(mRateArr):
        for y,maxRes in enumerate(sizeArr):
            if adaptiveTherapyBoolean:
                sequence = np.tile([False,True],int(days/2)+1)
                results[y,x]=tumorSizeCalcChangingMaxRes(sequence,mRate,numbOfSwitches,maxRes,n)
            else:
                sequence = np.full(days,True)
                results[y,x]=tumorSizeCalcChangingMaxRes(sequence, mRate,numbOfSwitches,maxRes,n)
    return results,mRateArr,sizeArr

switches = 10

#resAdaptive, xs, ys = varyMaxResStart(.01,.1,.01,1.001,2,switches,100,switches,True)
#resCont, xs, ys = varyMaxResStart(.01,.1,.01,1.001,2,switches,100,switches,False)

#resAdaptive, xs, ys = varyMRateSizeExp(.01,.1,.01,2,switches,100,True)
#resCont, xs, ys = varyMRateSizeExp(.01,.1,.01,2,switches,100,False)

resAdaptive, xs, ys = varyMRateSizeLinear(.01,.1,.01,2,switches,300,True)
resCont, xs, ys = varyMRateSizeLinear(.01,.1,.01,2,switches,300,False)


resDiff = resCont-resAdaptive
print(resDiff[0][0])

plt.imshow(resDiff)
maxVal = [-100,0,0]
for i in range(0,len(resDiff)):
    for j in range(0,len(resDiff[i])):
        if resDiff[i][j] > 0:
            print("switch = ", (i+2),"mRate =", (j+1)*.01,"difference = ", resDiff[i][j])
        if resDiff[i][j] > maxVal[0]:
            maxVal = [resDiff[i][j],i,j]
print(maxVal)
xAxis = np.arange(0, len(xs))
xs = [str(round(elem, 2)) for elem in xs]
xs = [elem.strip("0") for elem in xs]
plt.xticks(xAxis, xs)
plt.yticks(np.arange(0, len(ys)), ys)
plt.xlabel("Transition Rate")
plt.ylabel("Number of Switches")
plt.title("Tumor Size Difference (Continuous Therapy minus Pulsed Therapy)")
plt.colorbar()
plt.show()