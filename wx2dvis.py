import os
import wxdata as wxdata2
from pylab import *
from numpy import *

frames = 40
filename = "/home/sousae/dbg_warpx/recon_mmc/ssrecon_wv"

for n in range(0,frames+1):
		dh = wxdata2.WxData(filename, n)
		q = dh.read('qnew')
		x = linspace(q.grid.lowerBounds[0],q.grid.upperBounds[0],q.grid.numPhysCells[0])
		y = linspace(q.grid.lowerBounds[1],q.grid.upperBounds[1],q.grid.numPhysCells[1])
		X, Y = meshgrid(y,x)
		by = q[:, :, 0]
		
		pcolor(Y,X,by)
		
		draw()
