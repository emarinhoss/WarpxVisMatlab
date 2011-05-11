clear all; close all; clc;

%% Multilevel Monte Carlo
disp('Multilevel Monte Carlo')
L = 3;
Ly = 12.8;
frames = 40;


folder = '/home/sous776/scratch/mmc_recon/';
d = dir([folder 'recon_006*']);
n = size(d);

% allocate
count = zeros(1,3); s0=0; s1=0; s2=0;
flux0 = zeros(n(1),frames+1);
flux1 = zeros(n(1),frames+1);
flux2 = zeros(n(1),frames+1);
Ymean = zeros(L,frames+1);
Yvarn = Ymean;

disp('Processing directory:')
for k=1:n(1)
    
    disp(d(k).name)
    
    file0 = [folder d(k).name '/ssrecon_wv_0'];
    file1 = [folder d(k).name '/ssrecon_wv_1'];
    file2 = [folder d(k).name '/ssrecon_wv_2'];
    
    if (exist([folder d(k).name '/ssrecon_wv_2.pin'],'file') && ...
            exist([folder d(k).name '/ssrecon_wv_1.pin'],'file'))
        level = 2;
        s2 = s2 + 1;
    elseif (exist([folder d(k).name '/ssrecon_wv_1.pin'],'file'))
        level = 1;
        s1 = s1 + 1;
    else
        level = 0;
        s0 = s0 + 1;
    end
    
    switch level
        case 2   
            for l=0:frames
                [out0, x, y] = load_data_new(file0,'qnew',l);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux0(k,l+1) = dx*sum(abs(out0(mid,:,15)));
                
                [out1, x, y] = load_data_new(file1,'qnew',l);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux1(k,l+1) = dx*sum(abs(out1(mid,:,15)));
             
                [out2, x, y] = load_data_new(file2,'qnew',l);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux2(k,l+1) = dx*sum(abs(out2(mid,:,15)));
            end
            
            flux0(k,:) = flux0(k,:)/(2*Ly);
            flux0(k,:) = 0.2*flux0(k,:)/flux0(k,1);
            Ymean(1,:) = Ymean(1,:) + flux0(k,:);
            Yvarn(1,:) = Yvarn(1,:) + flux0(k,:).*flux0(k,:);

            flux1(k,:) = flux1(k,:)/(2*Ly);
            flux1(k,:) = 0.2*flux1(k,:)/flux1(k,1);
            Ymean(2,:) = Ymean(2,:) + flux1(k,:)-flux0(k,:);
            Yvarn(2,:) = Yvarn(2,:) + flux1(k,:).*flux1(k,:);
            
            flux2(k,:) = flux2(k,:)/(2*Ly);
            flux2(k,:) = 0.2*flux2(k,:)/flux2(k,1);            
            Ymean(3,:) = Ymean(3,:) + flux2(k,:)-flux1(k,:);
            Yvarn(3,:) = Yvarn(3,:) + flux2(k,:).*flux2(k,:);
            
        case 1
            for l=0:frames
                [out0, x, y] = load_data_new(file0,'qnew',l);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux0(k,l+1) = dx*sum(abs(out0(mid,:,15)));
                
                [out1, x, y] = load_data_new(file1,'qnew',l);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux1(k,l+1) = dx*sum(abs(out1(mid,:,15)));
            end
            
            flux0(k,:) = flux0(k,:)/(2*Ly);
            flux0(k,:) = 0.2*flux0(k,:)/flux0(k,1);
            Ymean(1,:) = Ymean(1,:) + flux0(k,:);
            Yvarn(1,:) = Yvarn(1,:) + flux0(k,:).*flux0(k,:);

            flux1(k,:) = flux1(k,:)/(2*Ly);
            flux1(k,:) = 0.2*flux1(k,:)/flux1(k,1);
            Ymean(2,:) = Ymean(2,:) + flux1(k,:)-flux0(k,:);
            Yvarn(2,:) = Yvarn(2,:) + flux1(k,:).*flux1(k,:);            
     
        otherwise
            for l=0:frames
                [out0, x, y] = load_data_new(file0,'qnew',l);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux0(k,l+1) = dx*sum(abs(out0(mid,:,15)));
            end
            
            flux0(k,:) = flux0(k,:)/(2*Ly);
            flux0(k,:) = 0.2*flux0(k,:)/flux0(k,1);
            Ymean(1,:) = Ymean(1,:) + flux0(k,:);
            Yvarn(1,:) = Yvarn(1,:) + flux0(k,:).*flux0(k,:);
    end
end

count(1) = s1+s2+s0;
count(2) = s1+s2;
count(3) = s2;

avg = zeros(1,frames+1); var=avg;
for m=1:L
    avg1 = Ymean(m,:)/count(m);
    var = var + 1/count(m)*(Yvarn(m,:)/count(m) - avg1.*avg1);
    
    avg = avg + avg1;
end
T = linspace(0,frames,frames+1);
errorbar(T,avg,sqrt(var),'r')