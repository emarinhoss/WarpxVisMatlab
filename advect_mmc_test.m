clear all; clc;

%% Multilevel Monte Carlo
L = 3;
frames = 40;


folder = '/home/sous776/WarpX/cfdlabruns/eder/PNNL/mmc_advect/';
d = dir([folder 'advect_001*']);
n = size(d);

% allocate
count = n(1)./[1 2 4];
Ymean = zeros(L,frames+1);
Yvarn = Ymean;

for k=1:n(1)
    file0 = [folder d(k).name '/advect_0'];
    for l=0:frames
        [o0,x,y] = advection_2d(file0,'qnew',l);
        maxP_0 = max(o0(round(length(y)/2),:));
        Ymean(1,l+1) = Ymean(1,l+1) + maxP_0;
        Yvarn(1,l+1) = Yvarn(1,l+1) + maxP_0*maxP_0;
    
        if exist([folder d(k).name '/advect_1.pin'])
            file1 = [folder d(k).name '/advect_1'];
            [o1,x,y] = advection_2d(file1,'qnew',l);
            maxP_1 = max(o1(round(length(y)/2),:));
            Ymean(2,l+1) = Ymean(2,l+1) + maxP_1-maxP_0;
            Yvarn(2,l+1) = Yvarn(2,l+1) + maxP_1*maxP_1;
        end
    
        if exist([folder d(k).name '/advect_2.pin'])
            file2 = [folder d(k).name '/advect_2'];
            [o2,x,y] = advection_2d(file2,'qnew',l);
            maxP_2 = max(o2(round(length(y)/2),:));
            %maxP_2-maxP_1
            Ymean(3,l+1) = Ymean(3,l+1) + maxP_2-maxP_1;
            Yvarn(3,l+1) = Yvarn(3,l+1) + maxP_2*maxP_2;
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
d2 = dir([folder 'advect_002*']);
n2 = size(d2);

for k=1:n2(1)
    file = [folder d2(k).name '/advect_mc'];
    for l=0:frames
        [o,x,y] = advection_2d(file,'qnew',l);
        maxP = max(o(round(length(y)/2),:));
        mc_avg(l+1) = mc_avg(l+1) + maxP;
        mc_var(l+1) = mc_var(l+1) + maxP*maxP;
    end
end

mc_avg = mc_avg/n2(1);
mc_var = mc_var/n2(1) - mc_avg.*mc_avg;

errorbar(0:frames,mc_avg,sqrt(mc_var),'b'), hold on
errorbar(0:frames,avg,sqrt(var),'r')