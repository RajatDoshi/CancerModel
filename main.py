from numba import njit
from scipy.integrate import odeint
import numpy as np
import gmpy2, math

class Model:
    def __init__(self,growthRates,mRatesUp,mRatesDown,resistanceArr,popArr=None,k=None,timesteps=None,drugWeeks=None):
        self.growthRates=growthRates
        self.mRatesUp=mRatesUp
        self.mRatesDown=mRatesDown
        self.resistances=resistanceArr
        self.popArr=popArr
        self.timesteps=timesteps
        self.drugWeeks=drugWeeks
        self.k=k

    def Run(self):
        change = np.zeros(self.popArr.size)
        return odeint(calculateChange, self.popArr, self.timesteps, (self.growthRates, self.mRatesUp, self.mRatesDown,
                                                                     self.k, self.resistances, self.drugWeeks, change),hmax=1)
    def RunEig(self, isDrug):
        return getEig(self.growthRates, self.resistances, self.mRatesUp, self.mRatesDown, isDrug)

def GenDefaultModel(mRate, nPops, step, days, drugSeq):
    popArr=calcPopDistribution(nPops,1)
    mRatesUp,mRatesDown=detMutationRate(mRate,nPops)
    return Model(np.full(nPops,.5), mRatesUp, mRatesDown, assignResistance(nPops,[0,1],0),
                 popArr, 1, np.arange(0, days, step), np.tile(drugSeq, int(days/(len(drugSeq)))+1))

def GenCustPopArr(mRate, popArr, step, days, drugSeq):
    nPops = len(popArr)
    mRatesUp,mRatesDown=detMutationRate(mRate,nPops)
    return Model(np.full(nPops,.25), mRatesUp, mRatesDown, assignResistance(nPops,[0,1],0),
                 popArr, 1, np.arange(0, days, step), periodicDrugDay(drugSeq,days))

def GenModelCustRes(mRate, nPops, resArr, odeStepSize, days, drugSeq):
    popArr = calcPopDistribution(nPops, 1)
    mRatesUp,mRatesDown=detMutationRate(mRate,nPops)
    return Model(np.full(nPops,.35), mRatesUp, mRatesDown, resArr,
                 popArr, 1, np.arange(0, days, odeStepSize), np.tile(drugSeq, int(days/len(drugSeq))+1))

def GenModelCustGrowthRate(mRate, nPops, resArr, growthRateArr, odeStepSize, days, drugSeq):
    popArr = calcPopDistribution(nPops, .1)
    mRatesUp,mRatesDown=detMutationRate(mRate,nPops)
    return Model(growthRateArr, mRatesUp, mRatesDown, resArr,
                 popArr, 1, np.arange(0, days, odeStepSize), np.tile(drugSeq, int(days/len(drugSeq))+1))

def calcPopDistribution(size,totalSize):
    popDistribution = np.zeros(size)
    for state in range(0,size):
        binomialCoef = gmpy2.bincoef(int(size-1),state)
        popDistribution[state] = totalSize*(binomialCoef/(math.pow(2,size)))
    if not  np.sum(popDistribution) == totalSize:
        factor = totalSize/np.sum(popDistribution)
        for i in range (0,len(popDistribution)):
            popDistribution[i] = popDistribution[i]*factor
    return popDistribution

#Purpose: Calculates Rate At Which the currState becomes finalState
#Variable: mRate = mutation rate
#Notes: Return Array Composed of all the rates
def transferRateArr(mRate,size):
    rateArrUp = np.zeros(size-1)
    rateArrDown = np.zeros(size-1)
    for currState in range(size-1):
        rateArrUp[currState]=(((size-1)-currState)/(size-1))*mRate
        rateArrDown[currState] = ((currState+1) / (size-1)) * mRate
    return rateArrUp,rateArrDown

#Purpose: Returns Mutation Rate Up and Down For Specific System Being Tested
#Variable: ourSystem and noTransfer are both booleans used to determine system
def detMutationRate(mRate,size,ourSystem=True,noTransfer=False):
    overallMRate = mRate
    if ourSystem:
        mutationRateUp,mutationRateDown = transferRateArr(overallMRate,size)
    elif noTransfer:
        mutationRateDown = [0]
        mutationRateUp = [0]
    elif size == 2:
        mutationRateUp = [overallMRate]
        mutationRateDown = [overallMRate]
    return mutationRateUp,mutationRateDown

def linearResistance(nPops, startPop, endPop, startRes, endRes):
    ret=np.zeros(nPops)
    if endPop-startPop>1:
        slope=(endRes-startRes)/((endPop-1) - startPop)
    for i in range(nPops):
        if(i<=startPop):
            ret[i]=startRes
        elif(i >= endPop):
            ret[i]=endRes
        else:
            ret[i]= slope * (i - startPop) + startRes
    return ret

#Purpose: Helper function to Calculate Linear Resistance
#Variable: size is the number of sub-populations, rangeArr is for resistance upper and lower bounds
def assignResistance(size,rangeArr,startPop=1):
    if size == 2:
        return [0,1]
    resistanceArr = np.zeros(size)
    lowerBound = rangeArr[0]; upperBound = rangeArr[1]
    step = (upperBound-lowerBound)/(size-startPop)
    for state in range(0,size):
        if state < startPop:
            resistanceArr[state] = 0
        else:
            resistanceArr[state] = lowerBound+((state-startPop+1) * step)
    return resistanceArr

#Purpose: Returns Resistance with exponential growth
#Variable: n is the growth factor constant that you can manipulate to accelerate resistance
def expResAssign(size,n,minIndex,maxIndex):
    resArr = np.zeros(size)
    for switch in range(0,size):
        if switch <= minIndex:
            resArr[switch] = 0
        elif switch >= maxIndex:
            resArr[switch] = 1
        else:
            resArr[switch] = (np.float128(n)**np.float128(switch-minIndex)-np.float128(1))/(np.float128(n)**np.float128(maxIndex-minIndex)-np.float128(1))
    return resArr

#Purpose: Incorporates a cost of resistances on growthRates
#Variable: maxDecreaseVal is how much percent lower should the maximum resistance population be
def costOfRes(growthRate,resistances,maxDecreaseVal=.2):
    growthRateArr = np.zeros(len(resistances))
    for index,res in enumerate(resistances):
        growthRateArr[index] = growthRate*(1-res*maxDecreaseVal)
    return growthRateArr

def costOfSwitches(nPops, minGR, maxGR):
    slope= (maxGR - minGR) / (nPops - 1)
    return [maxGR - slope * x for x in range(nPops)]

#Purpose: Changes drug sched everyday
#Variable: Sequence is what pattern of drug on and off you are going to repeat
def periodicDrugDay(sequence, days):
    return np.tile(sequence, int(days / len(sequence))+1)

def periodicDrugWeek(sequence, days):
    return np.tile(sequence, int(((days/7) / len(sequence))+1))

# Purpose: Given a population calculate mean and standard deviation
# Variable: popArr is population array composed of indivual subpopulation sizes in an array
# Note: mean and standard deviation are both weighted
def getMeanAndStandardDeviation(popArr):
    sizeOfPopArr = len(popArr)
    standardev = 0; numerator = 0
    for state in range(0, sizeOfPopArr):
        numerator = numerator + ((state / (sizeOfPopArr - 1)) * popArr[state])
    mean = numerator / sum(popArr)
    for state in range(0, sizeOfPopArr):
        standardev = standardev + ((((state / (sizeOfPopArr - 1)) - mean) ** 2) * popArr[state])
    standardev = standardev / sum(popArr)
    return mean, np.sqrt(standardev)

#Purpose: Determine Every Possible Treatment Schedule
#Variable: weeks is the total number of weeks administered
def therapyListGenerator(weeks):
    possibleCombinations = 2 ** weeks
    therapyArr = np.zeros((possibleCombinations,weeks+1))
    for i in range(0,possibleCombinations):
        for j in range(0,weeks+1):
            if j == weeks:
                therapyArr[i][j] = 0
            else:
                ans = int((i/(2**j))%2)
                if ans == 0:
                    therapyArr[i][j] = False
                else:
                    therapyArr[i][j] = True
    return therapyArr

#Purpose: Calculate Logistic Growth Rate for Cells
#Variables: gRate is growth rate, k is carrying capacity
@njit
def regrowthRate(gRate, subPopSize, totalSize, k):
    growth_rate = gRate*subPopSize*(1-(totalSize/k))
    return growth_rate

#Purpose: Calculate Death Rate due to Cancer Drug
#Variable: j is the Index Inside popArr Containing subPopulation
@njit
def deathRate(resistanceArr, drugPresent, index, subPopSize):
    if drugPresent > 0:
        deathRate = (1-resistanceArr[index]) * drugPresent * subPopSize
        return deathRate
    else:
        return 0

#Purpose:Calculate Change in Population Over Time for any System
#Variable: transitionRateUp/Down is absolute value of transition values in each direction
#Notes: for model with just sensitive and resistant --> [sensitive[drug,no drug], resistant[drug,no drug]]
@njit
def calculateChange(popArr,t,growthRates,mutationRateUp,mutationRateDown,k,resistanceArr,drugSched,change):
    totalPop=np.sum(popArr)
    for currState,subPopSize in enumerate(popArr):
        #calculating transition rates
        transitionRateFromUp=0
        transitionRateFromDown=0
        transitionRateUp=0
        transitionRateDown=0
        if currState<len(popArr)-1:
            transitionRateFromUp=mutationRateDown[currState]*popArr[currState+1]
            transitionRateUp = mutationRateUp[currState] * subPopSize
        if currState > 0:
            transitionRateFromDown = mutationRateUp[currState-1] * popArr[currState - 1]
            transitionRateDown = mutationRateDown[currState-1] * subPopSize
        #calculating growth and death rates
        birthTerm = regrowthRate(growthRates[currState], subPopSize, totalPop, k)
        if int(t)>len(drugSched)-1:
            drugOn=drugSched[len(drugSched)-1]
        else:
            drugOn=drugSched[int(t)]
        deathTerm = deathRate(resistanceArr, drugOn, currState, subPopSize)
        change[currState]=birthTerm - deathTerm - transitionRateUp - transitionRateDown +\
                          transitionRateFromUp + transitionRateFromDown
    return change

#Purpose: Calculates Eigen Values and Vector
#Variable: isDrug is a boolean to determine if drug is on or off
#Note: matrix,value,vector
def getEig(growthRates,resistanceArr,mutationRateUp,mutationRateDown,isDrug):
    nPops=len(growthRates);
    mat=np.zeros([nPops,nPops])
    for i in range(nPops):
        growthRate = 1 + growthRates[i]
        if i < nPops - 1:
          transitionRateFromUp = mutationRateDown[i]
          transitionRateUp = mutationRateUp[i]
          mat[i, i + 1] = transitionRateFromUp
          growthRate -= transitionRateUp
        if i > 0:
          transitionRateFromDown = mutationRateUp[i - 1]
          transitionRateDown = mutationRateDown[i - 1]
          mat[i,i-1]=transitionRateFromDown
          growthRate-=transitionRateDown
        if isDrug:
          growthRate-=(1-resistanceArr[i])
        mat[i,i]=growthRate
    eig=np.linalg.eig(mat); eigVals=eig[0]; vectors=eig[1];
    positives= np.min(vectors, axis=0) >= 0; negatives=np.max(vectors,axis=0) <= 0
    posVec=vectors[:,positives]; eigVal=eigVals[positives]
    if posVec.size==0:
        posVec=-vectors[:,negatives]
        eigVal=eigVals[negatives]
    eigVal=eigVal[0]
    posVec=np.reshape(posVec,[len(posVec)])
    return mat,eigVal,posVec