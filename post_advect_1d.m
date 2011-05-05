clear all; clc;

%% Multilevel Monte Carlo
disp('Multilevel Monte Carlo')
L = 3;
frames = 40;


folder = '/home/sous776/WarpX/cfdlabruns/eder/PNNL/mmc_advect/';
d = dir([folder 'advect_001*']);
n = size(d);

% allocate
count = n(1)./[1 2 4];
Ymean = zeros(L,frames+1);
Yvarn = Ymean;

disp('Processing directory:')
for k=1:n(1)
    
    disp(d(k).name)
    
    file0 = [folder d(k).name '/advect_0'];
    file1 = [folder d(k).name '/advect_1'];
    file2 = [folder d(k).name '/advect_2'];
    
    if (exist([folder d(k).name '/advect_2.pin'],'file') && ...
            exist([folder d(k).name '/advect_1.pin'],'file'))
        level = 2;
    elseif (exist([folder d(k).name '/advect_1.pin'],'file'))
        level = 1;
    else
        level = 0;
    end
    
    switch level
        case 2
            
            for l=0:frames
                [o0,x] = advection_1d(file0,l);
                maxP_0 = abs(x(2)-x(1))*sum(o0);
                Ymean(1,l+1) = Ymean(1,l+1) + maxP_0;
                Yvarn(1,l+1) = Yvarn(1,l+1) + maxP_0*maxP_0;
                
                [o1,x] = advection_1d(file1,l);
                maxP_1 = abs(x(2)-x(1))*sum(o1);
                Ymean(2,l+1) = Ymean(2,l+1) + maxP_1-maxP_0;
                Yvarn(2,l+1) = Yvarn(2,l+1) + maxP_1*maxP_1;
        
                [o2,x] = advection_1d(file2,l);
                maxP_2 = abs(x(2)-x(1))*sum(o2);
                Ymean(3,l+1) = Ymean(3,l+1) + maxP_2-maxP_1;
                Yvarn(3,l+1) = Yvarn(3,l+1) + maxP_2*maxP_2;
            end
            
        case 1
            
            for l=0:frames
                [o0,x] = advection_1d(file0,l);
                maxP_0 = abs(x(2)-x(1))*sum(o0);
                Ymean(1,l+1) = Ymean(1,l+1) + maxP_0;
                Yvarn(1,l+1) = Yvarn(1,l+1) + maxP_0*maxP_0;
                
                [o1,x] = advection_1d(file1,l);
                maxP_1 = abs(x(2)-x(1))*sum(o1);
                Ymean(2,l+1) = Ymean(2,l+1) + maxP_1-maxP_0;
                Yvarn(2,l+1) = Yvarn(2,l+1) + maxP_1*maxP_1;
            end
            
        otherwise
            for l=0:frames
                [o0,x] = advection_1d(file0,l);
                maxP_0 = abs(x(2)-x(1))*sum(o0);
                Ymean(1,l+1) = Ymean(1,l+1) + maxP_0;
                Yvarn(1,l+1) = Yvarn(1,l+1) + maxP_0*maxP_0;
            end
            
    end
end

avg = zeros(1,frames+1); var=avg; mc_avg=avg; mc_var=avg;
for m=1:L
    avg1 = Ymean(m,:)/count(m);
    var = var + 1/count(m)*(Yvarn(m,:)/count(m) - avg1.*avg1);
    
    avg = avg + avg1;
end

%% Monte Carlo
disp(' ')
disp('Monte Carlo')
d2 = dir([folder 'advect_002*']);
n2 = size(d2);

disp('Processing directory:')
for k=1:n2(1)
    file = [folder d2(k).name '/advect_mc'];
    disp(d2(k).name)
    for l=0:frames
        [o,x] = advection_1d(file,l);
        maxP = abs(x(2)-x(1))*sum(o);
        mc_avg(l+1) = mc_avg(l+1) + maxP;
        mc_var(l+1) = mc_var(l+1) + maxP*maxP;
    end
end

mc_avg = mc_avg/n2(1);
mc_var = mc_var/n2(1) - mc_avg.*mc_avg;

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

%% Exact calculation
disp(' ')
disp('Exact')
file = '/home/sous776/WarpX/cfdlabruns/eder/PNNL/mmc_advect/exact/advect_exact';

T = linspace(0,pi,frames+1); At=T;
for l=0:frames
    [o,x] = advection_1d(file,l);
    At(l+1) = abs(x(2)-x(1))*sum(o);
end

errorbar(T,mc_avg,sqrt(mc_var),'b'), hold on
errorbar(T,avg,sqrt(var),'r')
errorbar(T,pc_avg,pc_std,'k')
plot(T,At,'g')
legend('MC','MMC','PCM','Exact')
ylabel('Area')
xlabel('Time')