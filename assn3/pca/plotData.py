import numpy as np
import matplotlib.pyplot as plt 
import csv

data = []
with open('scatter.csv','rb') as f:
  points = csv.reader(f)
  for p in points:
    data.append([float(p[0]), float(p[1])])

data = np.array(data).transpose()
means = np.mean(data, axis=1)

for i in range(len(means)):
	data[i,:] -= means[i]
plt.axis([-1,1,-3,3])
plt.scatter(data[0,:],data[1,:])
plt.show()
