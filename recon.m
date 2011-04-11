clear all; clc;
component = 1;

filename='/home/sous776/WarpX/cfdlabruns/eder/recon/ssrecon_wv_freqcheck';
minval = 0;
maxval = 0;

for k=0:1
    [output, x, y] = load_data_new(filename,'qnew',k);
    pcolor(x,y,output(:,:,component)); colorbar;
    drawnow
    arrmin = min(min(output(:,:,component)));
    arrmax = max(max(output(:,:,component)));
    
    minval = min(minval,arrmin);
    maxval = max(maxval, arrmax);
end