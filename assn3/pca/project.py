import sys
import numpy as np
import matplotlib.pyplot as plt 
import csv

# if len(sys.argv) < 2:
#   print "Usage: python project.py <slope value>"
#   sys.exit(1)
# slope = float(sys.argv[1])

data = []
with open('scatter.csv','rb') as f:
	points = csv.reader(f)
 	for p in points:
 		data.append([float(p[0]), float(p[1])])

data = np.array(data).transpose()
means = np.mean(data, axis=1)


data[0,:] -= means[0]
data[1,:] -= means[1]

total_var = np.trace(np.cov(data))
curMax = 0
bestSlope = 0

for slope in np.arange(0.005, 6.0, .005):
	vec = np.array([0.0,1.0])
	if slope != np.inf:
	 	vec = np.array([1,slope])
	vec /= np.linalg.norm(vec)

	P = np.outer(vec,vec)
	#print "Projection Matrix\n",P

	slopeAligned = P.dot(data)

	phi = -np.arctan(slope)
	# clockwise rotation by phi
	R = np.array([[np.cos(phi), -np.sin(phi)],[np.sin(phi), np.cos(phi)]])
	xAligned = R.dot(slopeAligned)
	#compute variance here
	variance = np.var(xAligned[0,:])

	if variance > curMax:
		curMax = variance
		bestSlope = slope


print "Empirical Variance is %.4f"%(curMax)
print "Variance fraciton is %.4f"%(curMax/total_var)
print "Best slope is %.3f"%(bestSlope)
plt.scatter(slopeAligned[0,:],slopeAligned[1,:])
plt.show()

#subtract the first principle component from the dataset
#for i in range(len(vec)):
#	data[i,:] -= vec[i]

slope = -1.0/bestSlope

vec = np.array([1,slope])
vec /= np.linalg.norm(vec)

P = np.outer(vec,vec)
print "Projection Matrix\n",P

slopeAligned = P.dot(data)

phi = -np.arctan(slope)
# clockwise rotation by phi
R = np.array([[np.cos(phi), -np.sin(phi)],[np.sin(phi), np.cos(phi)]])
xAligned = R.dot(slopeAligned)
#compute variance here
variance = np.var(xAligned[0,:])


print "Empirical Variance is %.4f"%(variance)
print "Variance fraction is %.4f"%(variance/total_var)
print "Best slope is %.3f"%(slope)
plt.scatter(slopeAligned[0,:],slopeAligned[1,:])
plt.show()
