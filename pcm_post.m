frames = 40;
folder = '/home/sous776/WarpX/cfdlabruns/eder/PNNL/mmc_advect/';

%% Probabilist Collocation Method
disp(' ')
disp('PCM')
d2 = dir([folder 'advect_003*']);
n2 = size(d2);

pc_avg = zeros(1,frames+1); pc_std=pc_avg;

disp('Processing directory:')
for k=1:n2(1)
    file = [folder d2(k).name '/advect_pcm'];
    disp(d2(k).name)
    for l=0:frames
        fid = fopen([folder d2(k).name '/weight.w'],'r');
        wgt = fscanf(fid,'%f');
        fclose(fid);
        [o,x] = advection_1d(file,l);
        maxP = abs(x(2)-x(1))*sum(o);
        pc_avg(l+1) = pc_avg(l+1) + maxP*wgt;
        pc_std(l+1) = pc_std(l+1) + maxP*maxP*wgt;
    end
end

errorbar(0:40,pc_avg,pc_std)