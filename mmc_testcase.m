clear all; clc; close all;

%% Multilevel Monte Carlo
disp('Multilevel Monte Carlo')
L = 3;
frames = 10;
cells = 100;

folder = '/home/sous776/WarpX/cfdlabruns/eder/PNNL/mmc_advect/';
d = dir([folder 'advect_001*']);
n = size(d);

% allocate
count = n(1)./[1 2 4]; %zeros(1,L); s0=0; s1=1; s2=0;
Ymean = zeros(L,cells);
Yvarn = Ymean;

s = 2/pi;

disp('Processing directory:')
for k=1:n(1)
    
    disp(d(k).name)
    
    file0 = [folder d(k).name '/advect_0'];
    file1 = [folder d(k).name '/advect_1'];
    file2 = [folder d(k).name '/advect_2'];
    
    if (exist([folder d(k).name '/advect_2.pin'],'file') && ...
            exist([folder d(k).name '/advect_1.pin'],'file'))
        level = 2;
        %s2 = s2 + 1;
    elseif (exist([folder d(k).name '/advect_1.pin'],'file'))
        level = 1;
        %s1 = s1 + 1;
    else
        level = 0;
        %s0 = s0 + 1;
    end
    
    switch level
        case 2
            
            for l=frames:frames
                [o0,x0] = advection_1d(file0,l);
                maxP_0 = o0(2,:)./o0(1,:);
                Ymean(1,:) = Ymean(1,:) + maxP_0;
                Yvarn(1,:) = Yvarn(1,:) + maxP_0.*maxP_0;
                
                [o1,x1] = advection_1d(file1,l);
                maxP_1 = interp1(x1,o1(2,:)./o1(1,:),x0,'linear');
                Ymean(2,:) = Ymean(2,:) + maxP_1-maxP_0;
                Yvarn(2,:) = Yvarn(2,:) + maxP_1.*maxP_1;
        
                [o2,x2] = advection_1d(file2,l);
                maxP_2 = interp1(x2,o2(2,:)./o2(1,:),x0,'linear');
                Ymean(3,:) = Ymean(3,:) + maxP_2-maxP_1;
                Yvarn(3,:) = Yvarn(3,:) + maxP_2.*maxP_2;
            end
            
        case 1
            
            for l=frames:frames
                [o0,x0] = advection_1d(file0,l);
                maxP_0 = o0(2,:)./o0(1,:);
                Ymean(1,:) = Ymean(1,:) + maxP_0;
                Yvarn(1,:) = Yvarn(1,:) + maxP_0.*maxP_0;
                
                [o1,x1] = advection_1d(file1,l);
                maxP_1 = interp1(x1,o1(2,:)./o1(1,:),x0,'linear');
                Ymean(2,:) = Ymean(2,:) + maxP_1-maxP_0;
                Yvarn(2,:) = Yvarn(2,:) + maxP_1.*maxP_1;
            end
            
        otherwise
            for l=frames:frames
                [o0,x0] = advection_1d(file0,l);
                maxP_0 = o0(2,:)./o0(1,:);
                Ymean(1,:) = Ymean(1,:) + maxP_0;
                Yvarn(1,:) = Yvarn(1,:) + maxP_0.*maxP_0;
            end
            
    end
end

%count(1) = s0+s1+s2;
%count(2) = s1+s2;
%count(3) = s2;

avg = zeros(1,cells); var=avg; mc_avg=avg; mc_var=avg;
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
    for l=frames:frames
        [o,x] = advection_1d(file,l);
        maxP = o(2,:)./o(1,:);
        mc_avg = mc_avg + maxP;
        mc_var = mc_var + maxP.*maxP;
    end
end

mc_avg = mc_avg/n2(1);
mc_var = mc_var/n2(1) - mc_avg.*mc_avg;

%% Probabilist Collocation Method
disp(' ')
disp('PCM')
d2 = dir([folder 'advect_003*']);
n2 = size(d2);

pc_avg = zeros(1,cells); pc_var=pc_avg;

disp('Processing directory:')
for k=1:n2(1)
    file = [folder d2(k).name '/advect_pcm'];
    fid = fopen([folder d2(k).name '/weight.w'],'r');
    wgt = fscanf(fid,'%f');
    fclose(fid);
    disp(d2(k).name)
    for l=frames:frames
        [o,x] = advection_1d(file,l);
        maxP = o(2,:)./o(1,:);
        pc_avg = pc_avg + maxP*wgt;
        pc_var = pc_var + maxP.*maxP*wgt;
    end
end

%% Exact calculation
 disp(' ')
 disp('Exact')

T = linspace(0,1,cells); A=zeros(1,cells); B=A;

fid = fopen([folder 'mmc_values.txt'],'r');
mmc_vals = fscanf(fid,'%f');
fclose(fid); mmcval = mean(mmc_vals);

fid = fopen([folder 'mc_values.txt'],'r');
mc_vals = fscanf(fid,'%f');
fclose(fid);mcval = mean(mc_vals);

cs = sqrt(1.4);
wc = 10;
for l=0:9
    kn = 2*pi*(2*l+1);
    wn = sqrt(kn*kn*cs*cs+wc*wc);
    A = A - mmcval/(2*l+1)*sin(kn*T+wn*10.0);
    B = B - mcval/(2*l+1)*sin(kn*T+wn*10.0);
end

figure(1), p1=plot(T,mc_avg,'r',T,avg,'c',T,pc_avg,'--g',T,A,'--k');
legend('MC','MMC','PCM','Exact using MMC mean')
ylabel('u_x')
xlabel('Time')
title('MEAN')
set(p1,'LineWidth',2)


figure(2), p2=plot(T,mc_var,'r',T,var,'c',T,pc_var,'g');
legend('MC','MMC','PCM')
ylabel('\sigma^2')
xlabel('Time')
title('VARIANCE')
set(p2,'LineWidth',2)



figure(3),p3=plot(T,A-avg,'c',T,B-mc_avg,'r');
legend('MMC error','MC error')
ylabel('Error')
xlabel('X')
set(p3,'LineWidth',2)

%title('VARIANCE')