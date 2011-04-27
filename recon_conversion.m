clear all; clc; close all;

color = ['r','b','k','g','c','y','m','o','s'];

fileloc = '/home/sous776/WarpX/cfdlabruns/eder/PNNL/convergence/';
d = dir([fileloc 'res*']);
n = size(d);

frames = 40;
flux = zeros(n(1),frames+1);

for m=1:n(1)
    filename = [fileloc d(m).name '/ssrecon_wv_freqcheck'];
    
    for k=0:frames
        [output, x, y] = load_data_new(filename,'qnew',k);
        %subplot(2,1,1)
        %pcolor(x,y,output(:,:,component))
        %axis equal
        %colorbar
        %shading interp
        %drawnow
        dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
        flux(m,k+1) = dx * sum(abs(output(mid,:,15)));
    end
    Ly = 25.6;
    flux(m,:) = flux(m,:)/(2*Ly);
    flux(m,:) = 0.2*flux(m,:)/flux(m,1);
    plot(0:frames,flux(m,:),color(m)); hold on
    drawnow
end

title('Reconection Flux')
ylabel('Flux')
xlabel('Time')
legend(d(:).name)