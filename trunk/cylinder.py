import numpy
from pylab import *
from rbf_solver2 import rbf_solver

data = numpy.load('cylinder')

x = data['x']
y = data['y']
c = data['c']
r = data['r']
s = data['s']

rbf_solver(x, y, c, r, 0.007)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

x1   = arange(len(c))
  
ax1.scatter(x, y, c = c, marker = 'o', s = 10, edgecolors='none')    
ax2.plot(x1,c,'k--')
ax1.set_title('Solution')
ax2.set_title('point values')

plt.show()