clear all; clc;

filename='/home/sous776/WarpX/cfdlabruns/eder/PNNL/mmc_advect/advect_002_U_1.24729483555/advect_mc';

for i=1:20
    [output,x] = advection_1d(filename,i);
    plot(x,output)%,x,output(2,:),'r'), axis([-12.8 12.8 0 1])
    pause(0.1)
%    drawnow
end