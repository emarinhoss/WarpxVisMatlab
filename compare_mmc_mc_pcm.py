import glob
import os
import subprocess
import wxdata as wxdata2
from pylab import *
from numpy import *
import numpy


frames = 40
Ly = 12.8
wci = 0.1
T = linspace(0.0, frames*2.0*wci, frames+1)

########################################################################
########################################################################
# Multilevel Monte Carlo
########################################################################
########################################################################
mean_mmc = numpy.genfronmtxt("./mmc/mean.dat")
varn_mmc = numpy.genfronmtxt("./mmc/varn.dat")

########################################################################
########################################################################
# Monte Carlo
########################################################################
########################################################################

dmc = glob.glob("/home/sousae/scratch/mc/recon_001*")

flux0  = numpy.zeros((len(dmc),frames+1), numpy.float)
meanf0 = numpy.zeros(frames+1, numpy.float)
varf0  = numpy.zeros(frames+1, numpy.float)

for m in range(0,len(dmc)):
	filename = os.path.join(dmc[m], "ssrecon_wv")

	for n in range(0,frames+1):
		dh = wxdata2.WxData(filename, n)
		q = dh.read('qnew')
		dx = (q.grid.upperBounds[1]-q.grid.lowerBounds[1])/q.grid.numPhysCells[1]
		ny = int(round(q.grid.numPhysCells[1]/2))
		by = q[:, ny, 14]
		
		flux0[m,n] = dx*numpy.sum( numpy.fabs(by) )
		dh.close()

for l in range(0, size(d)):
	flux0[l,:] = 2*wci*flux[l,:]/flux0[l,0] # rescale to match GEM conditions
	meanf0 = meanf0 + flux0[l,:]
	varf0  = varf0 + flux0[l,:]*flux0[l,:]
	

varf0 = varf/(size(d)-1.0) - meanf*meanf/(size(d)*size(d)-size(d))
meanf0= meanf/size(d)

########################################################################
########################################################################
# Probabilistic Collocation Method
########################################################################
########################################################################

d = glob.glob("/home/sousae/scratch/pcm/recon_001*")

flux  = numpy.zeros((len(d),frames+1), numpy.float)
meanf = numpy.zeros(frames+1, numpy.float)
varf  = numpy.zeros(frames+1, numpy.float)
weight = []

for m in range(0,len(d)):
	filename = os.path.join(d[m], "ssrecon_wv")
	
	weigthfile=os.path.join(d[m], "weight.w")
	wts = open(weigthfile,'r')
	weight.append(float(wts.readline()))
	wts.close()
	
	for n in range(0,frames+1):
		dh = wxdata2.WxData(filename, n)
		q = dh.read('qnew')
		dx = (q.grid.upperBounds[1]-q.grid.lowerBounds[1])/q.grid.numPhysCells[1]
		ny = int(round(q.grid.numPhysCells[1]/2))
		by = q[:, ny, 14]
		
		flux[m,n] = dx*numpy.sum( numpy.fabs(by) )
		dh.close()
		
flux = flux/(2*Ly)

for l in range(0, size(d)):
	flux[l,:] = 2*wci*flux[l,:]/flux[l,0]* # rescale to match GEM conditions
	meanf = meanf + flux[l,:]*weight[l]
	varf  = varf + flux[l,:]*flux[l,:]*weight[l]

########################################################################
########################################################################
# Probabilistic Collocation Method
########################################################################
########################################################################

figure(1)
font = {'fontsize'   : 20}
subplot(2,1,1),plot(T,meanf,'g',T,meanf0,'r',mean_mmc[:,1],mean_mmc[:,1],'c',linewidth=2)
legend(('PCM','MC','MMC'),1)
xlabel('Time ($\omega_{ci}t$)',font)
ylabel('Reconnected Flux',font)
yticks(fontsize=18)
xticks(fontsize=18)

subplot(2,1,2),plot(T,varnf,'g',T,varnf0,'r',varn_mmc[:,1],varn_mmc[:,1],'c',linewidth=2)
legend(('PCM','MC','MMC'),1)
xlabel('Time ($\omega_{ci}t$)',font)
ylabel('Variation, $\sigma^2$',font)
yticks(fontsize=18)
xticks(fontsize=18)

savefig('recon_comparison.png')
