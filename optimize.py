from scipy.integrate import odeint
import numpy as np
import main
from scipy.optimize import minimize
import random

#Purpose: Create Objective Function
#Variable: sequence is the schedule that the medication is administered
def objectiveFun(sequence):
    days = 10; size = 5; k = 1; mRate = .05; growthRates = np.full(size,.5)
    popArr = main.calcPopDistribution(size,.01); resistances = main.expResAssign(size,2,0,size-1)
    mutationRatesUp, mutationRatesDown = main.detMutationRate(mRate,size,True,False)
    resultArr = odeint(main.calculateChange, popArr, np.arange(0, days, 1), (growthRates, mutationRatesUp, mutationRatesDown, k, resistances, main.periodicDrugDay(sequence, days),np.zeros(popArr.size)))
    weightAns = 0
    for i in range(0, len(resistances)):
        for j in range(0,len(resultArr)):
            weightAns = weightAns + (resistances[i]*resultArr[j][i])
    return weightAns

def constraintOne(sequence):
    return sum(sequence)-1

size = 10; aguess = np.full(int(size/2),1); bguess = np.full(int(size/2),0); initGuess = np.append(aguess,bguess)
b = (0,1); bnds = (b,b)*(int(size/2))
con1 = {'type': 'ineq', 'fun': constraintOne}
solution = minimize(objectiveFun,initGuess,method='SLSQP',bounds=bnds,constraints=con1)

while solution.fun > .005:
    print(solution.fun); print(solution.x); print("")
    for index in range(0,size):
        initGuess[index] = random.randint(0,1)
    b = (0, 1); bnds = (b, b) * (int(size / 2))
    con1 = {'type': 'ineq', 'fun': constraintOne}
    solution = minimize(objectiveFun, initGuess, method='SLSQP', bounds=bnds, constraints=con1)

oneCounter = 0; zeroCounter = 0
therapyComb = np.array(solution.x)
for i in range(0, len(therapyComb)):
    if therapyComb[i] > 0:
        therapyComb[i] = 1
        oneCounter = oneCounter + 1
    else:
        therapyComb[i] = 0
        zeroCounter = zeroCounter + 1
print("Zero", zeroCounter); print("One:", oneCounter)
print(therapyComb)
print(solution.fun)