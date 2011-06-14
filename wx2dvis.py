import os
import subprocess
import wxdata as wxdata2
from pylab import *
from numpy import *


frames = 40
location = "/home/sousae/dbg_warpx/recon_mmc/ssrecon_wv"

for n in range(0,frames+1):
		dh = wxdata2.WxData(location, n)
		q = dh.read('qnew')
		x = linspace(q.grid.lowerBounds[0],q.grid.upperBounds[0],q.grid.numPhysCells[0])
		y = linspace(q.grid.lowerBounds[1],q.grid.upperBounds[1],q.grid.numPhysCells[1])
		X, Y = meshgrid(y,x)
		by = q[:, :, 0]
		
		pcolor(Y,X,by)
		filename = str('%03d' % n) + '_recon.png'
		savefig(filename, dpi=100)
		clf()
		
		dh.close()


command = ('mencoder',
           'mf://*_recon.png',
           '-mf',
           'type=png:w=800:h=600:fps=24',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')
           
subprocess.check_call(command)           
