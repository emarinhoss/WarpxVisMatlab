clear all; clc;

filename='/home/sous776/WarpX/cfdlabruns/eder/recon/ssrecon_wv_freqcheck';
frames = 40;
component = 1;

flux = zeros(1,frames+1);

for k=0:frames
    [output, x, y] = load_data_new(filename,'qnew',k);
    subplot(2,1,1)
    pcolor(x,y,output(:,:,component))
    %axis equal
    colorbar
    shading interp
    drawnow
    
    dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
    flux(k+1) = dx * sum(abs(output(mid,:,15)));
end
Ly = y(end)-y(1);
flux = flux/(2*Ly);
flux = 0.2*flux/flux(1);

subplot(2,1,2),plot(0:40,flux);
title('Reconection Flux');
ylabel('Flux');
xlabel('Time');