clear all; clc;

filename='/home/sousae/WarpX/cfdlabruns/eder/advection/testcase/advection_1d_test';

for i=1:100
    [output,x] = advection_1d(filename,'qnew',i);
    [reint,xi] = advection_1d(filename,'qnint',i);
    subplot(3,1,1),plot(x,output(1,:),'b','LineWidth',2), axis([-12.8 12.8 0 1])
    title('Interior/Computational Domain')
    subplot(3,1,2),plot(x,output(2,:),'r','LineWidth',2), axis([-12.8 12.8 0 1])
    title('Artificial Boundary Condition')
    subplot(3,1,3),plot(xi,reint,'k','LineWidth',2), axis([-12.8 12.8 0 1])
    xlabel('Position')
    title('ABC Re-integrated Value')
    pause
%    drawnow
end