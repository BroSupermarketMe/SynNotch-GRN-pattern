% screen_params: 参数扫描：在自然生长条件下（check_params），
% 使用不同参数，获取稳态后的各项指标（Dmax/Dmin，达到稳态的时间，Gal4和GFP的H值
% （非活性细胞周围有多少个活性细胞）），并保存结果到MatData/result.mat以及MatData/params.mat

clear all;
clc

% rng(20240110);
params.rho = 30;

params.tspan=[0 5];
params.Cellcyclelength = 1;
params.max_Cellcyclelength = 2;
params.N = 8; % N>0

params.beta_N = logspace(-1, 3, params.rho);
params.beta_G = logspace(-1, 3, params.rho);
params.beta_F = logspace(-1, 3, params.rho);

% params.beta_N= 10^(0.9655);
% params.beta_G= 10^(1.2759);
% params.beta_F= 10^(0.0345);

params.mu = 1;
params.nu = 1;
params.xi = 1;
params.k1 = 1;

params.n1 = 4;
params.n2 = 4;

entropy = zeros(params.rho,params.rho,params.rho);
surround_surEGFP = zeros(params.rho,params.rho,params.rho);
surround_Gal4 = zeros(params.rho,params.rho,params.rho);
Dmax = zeros(params.rho,params.rho,params.rho);
Dmin = zeros(params.rho,params.rho,params.rho);
T = zeros(params.rho,params.rho,params.rho);

iterations = size(entropy,[1,2,3]);

ppm = ParforProgressbar(prod(iterations));
parfor ix = 1:prod(iterations)
    [para_i,para_j,para_k] = ind2sub(iterations,ix);
    paras = params;
    paras.beta_N = params.beta_N(para_i);
    paras.beta_G = params.beta_G(para_j);
    paras.beta_F = params.beta_F(para_k);
    [yout,tout, entropy(ix), surround_surEGFP(ix), surround_Gal4(ix)] = check_params(paras);
    ppm.increment();

    % finding max and min values of D
    k=size(yout,2)/4;
    Dmax(ix)=max(max(yout(end,2*k+1:3*k)));
    Dmin(ix)=abs(min(min(yout(end,2*k+1:3*k))));
    
    % finding cases where patterning occurs (when Dmax/Dmin>1.2)
    % and getting the patterning time 
    if Dmax(ix)/Dmin(ix)>1.2 
        T(ix)=getPatterningTime(tout,yout,k,Dmax(ix),Dmin(ix));
    else
        T(ix)=NaN; % patterning time is not set for the no patterning case
    end
end
delete(ppm);

save('MatData/result.mat','entropy','surround_surEGFP','surround_Gal4','T','Dmax','Dmin')
save('MatData/params.mat','params')


function T=getPatterningTime(tout,yout,k,Dmax,Dmin)
    
    % This function estimates patterning time. This is done by the follwoing 3 
    % steps:
    % 1. find all the high D cells ('onCells')
    % 2. find the time it takes for each 'on cell' to reach 90% of its final
    % level ('TonCells')
    % 3. get median value of the times calculated in stage 2
    
    onCells=find(yout(end,2*k+1:3*k)>0.5*(Dmax+Dmin));
    for i=1:length(onCells)
        Tind=find(yout(:,2*k+onCells(i))>0.9*yout(end,2*k+onCells(i)),1,'first');
        TonCells(i)=tout(Tind);
    end 
    
    T=median(TonCells);

end