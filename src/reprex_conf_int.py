import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

import scikits.bootstrap as bootstrap

lowess = sm.nonparametric.lowess

expression = np.random.uniform(3, 5,size=(100))
pseudotime = np.random.uniform(0, 1,size=(100))

fig, ax = plt.subplots()

plt.scatter(pseudotime, expression)

Replications = np.array([np.random.choice(expression, size=len(expression), replace=True) for _ in range(1000)])
Mean = np.mean(Replications, axis=1)

lowb=abs(Mean-np.percentile(Mean,5,interpolation='nearest')).argmin()  
upb=abs(Mean-np.percentile(Mean,95,interpolation='nearest')).argmin() 

lowc=Replications[lowb,:]
upc=Replications[upb,:]

z = lowess(expression, pseudotime)
zlow = lowess(lowc, pseudotime)
zhi = lowess(upc, pseudotime)

ax.plot(z[:,0], z[:,1], color="gray")
ax.plot(zlow[:,0], zlow[:,1], color="red")
ax.plot(zhi[:,0], zhi[:,1], color="red")

plt.show()
