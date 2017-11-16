import matplotlib.pyplot as plt
from math import pi
from numpy import arange, cos, degrees

## TODO:
#  Play around with the weights defined in k
#  and the offsets defined in offset. See
#  how the form of the potential energy, U, changes.

#k = {1:1, 2:5.0, 3:1}
#offset = {1:-1*pi, 2:0, 3:pi}

k = {1:2.0, 2:1.0, 3:2.0}
offset = {1:-1*pi, 2:0, 3:pi}

#k = {1:2.0, 2:1.0, 3:2.0}
#offset = {1:-1*pi, 2:0, 3:pi}



x = arange(0,2*pi,0.05)
U = sum((k[n]*(1-cos(n*x-offset[n])) for n in (1,2,3)))

plt.plot(degrees(x),U)
plt.xlim(0, 360)
plt.ylabel('Potential Energy')
plt.xlabel('Angle')
plt.xticks(range(0, 361, 60))
plt.show()
