import numpy as np 
import matplotlib.pyplot as plt
import sys 

#name="file"+sys.argv[1]
name="vector0-79"
data=np.loadtxt(name)
plt.hist(data)
plt.show()
