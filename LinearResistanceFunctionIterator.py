from main import Model
import numpy as np
import os

maxSwitches=11
minSwitches=1
ResSubDivs=10
contTherapy=[True]
pulsedTherapy=[True,False]
outputPath="~"
fullPath=outputPath+"LinearResistancePlots"
os.makedirs(fullPath)

for startX in range(minSwitches,maxSwitches):
  for startY in range(ResSubDivs-1):
    nCases=0
    for endX in range(startX+1,maxSwitches+1):
      for endY in range(startY+1,ResSubDivs):
        nCases+=1
    print(nCases)
