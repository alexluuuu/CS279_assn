# You are to implement the stochastic and laplacian diffusers in the code
# below. Specifically, you are to to implement the update(self) functions in
# the StochasticDiffuser and LaplacianDiffuser classes. Once you implement the
# update functions, you can invoke the script as follows:
#
# $ python diffuser.py stoc OR $ python diffuser.py lap
#
# Note use 'stoc' for the stochastic diffuser and 'lap' for the laplacian
# diffuser. Change the initial settings and parameters for the simulation
# below in the main method. The units in the StochasticDiffuser are number of
# particles. The units in the LaplacianDiffuser are an aribitrary measure of
# concentration.

from numpy import zeros, max as npmax
from random import random, randint
from math import floor,ceil

class Diffuser(object):

  def __init__(self,
               nX = 25,
               nY = 25,
               diffusion = 1,
               ):
    self.nX = nX
    self.nY = nY
    self.blocks = {(x,y):0 for x in xrange(nX) for y in xrange(nY)}
    self.diffusion = diffusion

  def setBlock(self,block,val):
    if block in self.blocks:
      self.blocks[block] = val

  def getBlock(self,block):
    if block in self.blocks:
      return self.blocks[block]
    return None

  def getGrid(self):
    grid = zeros((self.nX,self.nY))
    for (x,y), vals in self.blocks.iteritems():
      grid[x,y] = vals
    return grid

  def simulateTimeStep(self,nIters=1):
    for step in xrange(nIters):
      self.blocks = self.update()
      

  def update(self):
    raise NotImplementedError("Only instantiate Stochastic or Laplacian Diffusers")


  # Returns the correct x index after moving step blocks
  # to the right (to the left for negative step)
  def nextX(self,x,step):
    return (x+step)%self.nX
  # Returns the correct y index after moving step blocks
  # up (down for negative step)
  def nextY(self,y,step):
    return (y+step)%self.nY


class StochasticDiffuser(Diffuser):

  def update(self):
    # First, all values in new blocks are set to 0
    newBlocks = {xy:0 for xy in self.blocks} 

    # Next we will iterate through the keys (i.e. the coordinates) and values
    # (containing the number of particles) of self.blocks
    for (i,j), block in self.blocks.iteritems(): 

      # The delta/step size depend on the diffusion constant which is set in
      # the main method
      delta = ceil(self.diffusion)

      # Now we will iterate through every particle in each block
      for part in range(block): 
        # Introduce some randomness in step size
        step = delta - (random() > 0.7)
        direction = randint(1,4)

        if direction == 1:
          newBlocks[(self.nextX(i, step), j)] += 1
        elif direction == 2:
          newBlocks[(self.nextX(i, -1*step), j)] += 1

        elif direction == 3:
          newBlocks[(i, self.nextY(j, step))] += 1
        elif direction == 4:
          newBlocks[(i, self.nextY(j, -1*step))] += 1
    
    return newBlocks


class LaplacianDiffuser(Diffuser):


  def ThreePointApprox(self, curX, curY, grid):
    coefficients = [1.0 , -2.0, 1.0]
    dx2 = 0
    dy2 = 0
    h = 1.0
    for i, c in enumerate(coefficients, -1):
      dx2 += c*grid[(self.nextX(curX, i*h), curY)]
      dy2 += c*grid[(curX, self.nextY(curY, i*h))]

    return (dx2, dy2)

  def update(self):
    # First, all values in the derivative matrices are set to 0
    ddx = {xy:0 for xy in self.blocks}
    ddy = {xy:0 for xy in self.blocks}
    ddt = {xy:0 for xy in self.blocks}
    temp_blocks = self.blocks
    h = 1.0
    # Next we will iterate through the keys of self.blocks (i.e. the
    # coordinates) to compute the first discrete derivatives
    for (i,j) in self.blocks: 
     
      ddx[(i,j)] = (1.0/(h))*( self.blocks[(self.nextX(i, h), j)] - self.blocks[(self.nextX(i, 0), j)] )
      ddy[(i, j)] = (1.0/(h))*( self.blocks[(i, self.nextY(j, h))] - self.blocks[(i, self.nextY(j, 0))] )
    
    # Finally, we will once again iterate through the keys of self.blocks to
    # compute the second discrete derivatives
    for (i,j) in self.blocks:
      
      lapx, lapy = self.ThreePointApprox(i, j, self.blocks)
      ddt[(i,j)] = self.diffusion*(lapx + lapy)
    
    for x in ddt:
      temp_blocks[x] = self.blocks[x] + ddt[x]
    return temp_blocks



def main():
  from matplotlib.pyplot import *
  from matplotlib.animation import FuncAnimation
  from numpy import array
  from matplotlib import colors
  import sys

  usage = "Usage: python diffuser.py [stoc | lap]"
  iters = 10000
  cm = "cool"
  length = 25

  if len(sys.argv) < 2:
    print usage
    exit(1)


  # EDIT DIFFUSION PARAMETERS below for Stochastic and Laplacian Diffusers
  if sys.argv[1] == "stoc":
    # play around with different values for diffusion
    # (it only makes sense to use integers here)
    diffusion = 2
    c = StochasticDiffuser(nX = length,
                           nY = length,
                           diffusion = diffusion,
                           )
  elif sys.argv[1] == "lap":
    # play around with different values for diffusion
    # (note what happens when diffusion > 0.25)
    diffusion = .24
    c = LaplacianDiffuser(nX = length,
                          nY = length,
                          diffusion = diffusion,
                          )
  else:
    print usage
    exit(1)


  """ Here are a variety of settings for the initial conditions of the
      diffusion simulations.  Feel free to try some of your own. Note:  The
      units in the StochasticDiffuser are number of particles.  The units in
      the LaplacianDiffuser are an aribitrary measure of Concentration.  """

  ##############################################################
  # EDIT INITIAL SIMULATION PARAMETERS BELOW. 

  """SINGLE PARTICLE"""
  #c.setBlock((length/2,length/2),1)

  """Point Mass"""
  c.setBlock((length/2,length/2),1250)

  """1D Diffusion"""
  # for i in range(length):
  #   c.setBlock((i,length/2),125)

  """1D Gradient Diffusion"""
  # for i in range(length):
  #   for j in range(length):
  #     c.setBlock((i,j), 2*abs(length/2-i))

  # END PARAMTERS
  ##############################################################


  # This portion of the code runs the simulation through matplotlib
  prev_grid = zeros((length, length))
  for l in range(iters):
    clf()
    grid = c.getGrid()
    #if abs(npmax(grid) - npmax(prev_grid)) < .001:
    #  print l
    #  print 'converged!'
    #  break
    pcolor(array([i for i in xrange(length+1)]),array([j for j in xrange(length+1)]),grid,cmap = cm)
    colorbar()
    #un-comment the line below to 'fix' the colorbar to the range [0,3] so it is not longer dynamic
    #clim(0,3) 
    xlabel("Iteration %d"%l)
    pause(.01)
    c.simulateTimeStep(1)
    prev_grid = grid


# Boiler plate invokes the main function
if __name__ == "__main__":
  main()

