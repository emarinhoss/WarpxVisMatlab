clear all; clc; close all;

%color = ['r','b','k','g','c','y','m','o','s'];

fileloc = '/home/sous776/WarpX/cfdlabruns/eder/PNNL/NERCS_RUNS/';
d = dir([fileloc 'recon_005*']);
n = size(d);

frames = 40;
flux = zeros(n(1),frames+1);

for m=1:n(1)
    filename = [fileloc d(m).name '/ssrecon_wv'];
    
    for k=0:frames
        [output, x, y] = load_data_new(filename,'qnew',k);
        dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
        flux(m,k+1) = dx * sum(abs(output(mid,:,15)));
    end
    Ly = 25.6;
    flux(m,:) = flux(m,:)/(2*Ly);
    flux(m,:) = 0.2*flux(m,:)/flux(m,1);
end

meanf = sum(flux)/n(1);
varf = sum(flux.*flux)/n(1) - meanf;
stdf = sqrt(abs(varf));

errorbar(0:40,meanf,varf);
title('Mean and variance of the reconnection flux for 16 MC runs.')
ylabel('Flux')
xlabel('Time')
%legend(d(:).name)