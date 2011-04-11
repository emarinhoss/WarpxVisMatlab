clear all; clc;

filename='/home/sousae/WarpX/cfdlabruns/eder/advection/testcase/advection_1d_test';

for i=1:100
    [output,x] = advection_1d(filename,i);
    plot(x,output(1,:)), axis([-12.8 12.8 0 0.004])
    pause(0.2)
%    drawnow
end