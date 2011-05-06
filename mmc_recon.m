%% Multilevel Monte Carlo
disp('Multilevel Monte Carlo')
L = 3;
frames = 40;


folder = '/home/sous776/scratch/mmc_recon/';
d = dir([folder 'recon_006*']);
n = size(d);

% allocate
s = 0.2/(2*12.8);   %sclaling
count = n(1)./[1 2 4];
flux0 = zeros(count(1),frames+1);
flux1 = zeros(count(2),frames+1);
flux2 = zeros(count(3),frames+1);
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
    elseif (exist([folder d(k).name '/ssrecon_wv_1.pin'],'file'))
        level = 1;
    else
        level = 0;
    end
    
    switch level
        case 2
            
            for l=0:frames
                [output, x, y] = load_data_new(file0,'qnew',k);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux0(k,l+1) = s*dx * sum(abs(output(mid,:,15)));
                Ymean(1,l+1) = Ymean(1,l+1) + flux0(k,l+1);
                Yvarn(1,l+1) = Yvarn(1,l+1) + flux0(k,l+1)*flux0(k,l+1);
                
                [output, x, y] = load_data_new(file0,'qnew',k);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux1(k,l+1) = s*dx * sum(abs(output(mid,:,15)));
                Ymean(2,l+1) = Ymean(2,l+1) + flux1(k,l+1)-flux0(k,l+1);
                Yvarn(2,l+1) = Yvarn(2,l+1) + flux1(k,l+1)*flux1(k,l+1);
                
                [output, x, y] = load_data_new(file0,'qnew',k);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux2(k,l+1) = s*dx * sum(abs(output(mid,:,15)));
                Ymean(3,l+1) = Ymean(3,l+1) + flux2(k,l+1)-flux1(k,l+1);
                Yvarn(3,l+1) = Yvarn(3,l+1) + flux2(k,l+1)*flux2(k,l+1);
            end
            
        case 1
            
            for l=0:frames
                [output, x, y] = load_data_new(file0,'qnew',k);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux0(k,l+1) = s*dx * sum(abs(output(mid,:,15)));
                Ymean(1,l+1) = Ymean(1,l+1) + flux0(k,l+1);
                Yvarn(1,l+1) = Yvarn(1,l+1) + flux0(k,l+1)*flux0(k,l+1);
                
                [output, x, y] = load_data_new(file0,'qnew',k);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux1(k,l+1) = s*dx * sum(abs(output(mid,:,15)));
                Ymean(2,l+1) = Ymean(2,l+1) + flux1(k,l+1)-flux0(k,l+1);
                Yvarn(2,l+1) = Yvarn(2,l+1) + flux1(k,l+1)*flux1(k,l+1); 

            end
            
        otherwise
            for l=0:frames
                [output, x, y] = load_data_new(file0,'qnew',k);
                dx = abs(x(2)-x(1)); mid=ceil(0.5*length(y));
                flux0(k,l+1) = s*dx * sum(abs(output(mid,:,15)));
                Ymean(1,l+1) = Ymean(1,l+1) + flux0(k,l+1);
                Yvarn(1,l+1) = Yvarn(1,l+1) + flux0(k,l+1)*flux0(k,l+1);
            end
            
    end
end

avg = zeros(1,frames+1); var=avg;
for m=1:L
    avg1 = Ymean(m,:)/count(m);
    var = var + 1/count(m)*(Yvarn(m,:)/count(m) - avg1.*avg1);
    
    avg = avg + avg1;
end

errorbar(T,avg,sqrt(var),'r')