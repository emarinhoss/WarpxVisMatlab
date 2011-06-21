import os
import subprocess
import wxdata as wxdata2
from pylab import *
import numpy as np
from numpy import *


frames = 40
location = "/home/sousae/dbg_warpx/recon_mmc/ssrecon_wv"
minpress = 10.0

for n in range(0,frames+1):
		dh = wxdata2.WxData(location, n)
		q = dh.read('qnew')
		x = np.linspace(q.grid.lowerBounds[0],q.grid.upperBounds[0],q.grid.numPhysCells[0])
		y = np.linspace(q.grid.lowerBounds[1],q.grid.upperBounds[1],q.grid.numPhysCells[1])
		X, Y = np.meshgrid(y,x)
		by = q[:, :, 4]
		
		minpress = min(minpress,by.min())
		pcolor(Y,X,by), colorbar()
		filename = str('%03d' % n) + '_recon.png'
		savefig(filename, dpi=100)
		clf()
		
		dh.close()


command = ('mencoder',
           'mf://*_recon.png',
           '-mf',
           'type=png:w=800:h=600:fps=10',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')
           
subprocess.check_call(command)
f = open('pressure_minimum', 'w')
f.write(str(minpress))
f.close()
os.system('rm *_recon.png')
