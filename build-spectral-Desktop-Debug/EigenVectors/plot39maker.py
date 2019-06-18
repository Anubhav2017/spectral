import numpy as np 
import matplotlib.pyplot as plt
import sys 

name="file39"
data=np.loadtxt(name)

N=data.shape[0]

mindata=-0.039
maxdata=0.028
datarange=maxdata-mindata
relativevalues=[]
for i in range(N):
    value=0
    if(data[i]<mindata):
        value=0
    elif(data[i]>maxdata):
        value=1
    else:
        value=(data[i]-mindata)/datarange
    relativevalues.append(value)

plt.hist(relativevalues)
plt.show()
