clear all; clc; close all;

fileloc = '/home/sousae/scratch/mmc_recon/';
d = dir([fileloc 'recon_007*']);
n = size(d); num = n(1);

frames = 40;
flux = zeros(num,frames+1);

for m=1:num
    clc
    filename = [fileloc d(m).name '/ssrecon_wv']
    m
    for k=0:frames
        k
        [output, x, y] = load_data_new(filename,'qnew',k);
        dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
        flux(m,k+1) = dx * sum(abs(output(mid,:,15)));
    end
    Ly = 12.8;
    flux(m,:) = flux(m,:)/(2*Ly);
    flux(m,:) = 0.2*flux(m,:)/flux(m,1);
end

meanf = sum(flux)/num;
varf = sum(flux.*flux)/num - meanf.*meanf;
stdf = sqrt(abs(varf));

errorbar(0.1*(0:40),meanf,varf), hold on;
%plot(0:40,flux(1,:),'--g',0:40,flux(num,:),'--r')
title('Magnetic Reconnection Flux ')
ylabel('Flux')
xlabel('Time')
%legend(d(:).name)
