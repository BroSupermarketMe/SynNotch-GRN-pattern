% Copyright (C) 2014, David Sprinzak
% This program is part of Lateral Inhibition Tutorial.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% multicell_sysNotch：周期性边界条件下的多细胞模型。

function [yout,tout,H,params,F] = multicell_sysNotch(params, y0, ifplot)

% multicell_LI simulates lateral inhibition in a hexagonal lattice. 
% The structure params contains the model parameters of the system. 
% TOUT is a vector containing the time points of the solution 
% between 0 and Tmax. YOUT is a matrix containing the numerical 
% solution for each variable for each time point. Each row in 
% YOUT is a vector of the size of TOUT. F is a movie of the simulation. 
% get the default parameters if none provided
arguments
    params struct = defaultparams
    % setting the initial conditions + noise
    y0 (:,1) double = getIC(params, params.P*params.Q)
    ifplot logical = true
end

% set the random number generator
% rng(1);

Tmax=90; tspan=[0 Tmax]; % set time for simulation

P=params.P;     % number of cells per column
Q=params.Q;     % number of columns - MUST BE EVEN
k=P*Q;          % number of cells

% get the connectivity matrix
params.connectivity=getconnectivityM(P,Q);

% run simulation with lateral inhibition
options = odeset('nonnegative',1:4*k);  % make sure all concentrations are non-negative
[tout,yout] = ode45(@(t,y) li(t,y,params),tspan,y0,options);

% get the criterion for lateral inhibition
H1=getCriterion(yout(end,2*k+1:3*k),P,Q,k,32);
H2=getCriterion(yout(end,3*k+1:4*k),P,Q,k,32);
H=[H1,H2];

if ifplot
    % show time traces of two cells with lateral inhibition
    plot2cells(tout,yout,k);
    % 
    % show lattice simulation
    F=movielattice(tout,yout(:,2*k+1:3*k),P,Q,k);
    
    % % Hexagonal FFT for surEGFP
    % addpath('HFFT/');
    % X=HexImage_pq(yout(end,2*k+1:3*k),P,Q);
    % padsize=4;
    % X=padding(X,padsize/2);
    % yout_padding = PQImage_hex(X);
    % movielattice(tout(end,:),abs(yout_padding),P+padsize,Q+padsize,(P+padsize)*(Q+padsize));
    % F=hfft2(X);
    % F=hfftshift(F);
    % yout_f = PQImage_hex(F);
    % % 绘制横轴为空间频率的相位，纵轴为幅度的图像，只绘制非零幅度的实心点
    % F_n=get_X_n(F);
    % figure,plot(mod(rad2deg(angle(F_n(:))),360),abs(yout_f),'.'),xlim([1,360]);
    % xlabel('旋转角（度）')
    % ylabel('幅度')
    % movielattice(tout(end,:),abs(yout_f),P+padsize,Q+padsize,(P+padsize)*(Q+padsize));

    % save("./result/fft_result.mat", 'yout_f','tout','-mat');
    % yout_f = remove_padding("./result/fft_result.mat", "./result/fft_background.mat");
    % figure,plot(mod(rad2deg(angle(F_n(:))),360),abs(yout_f),'.'),xlim([1,360]);
    % movielattice(tout(end,:),abs(yout_f),P+padsize,Q+padsize,(P+padsize)*(Q+padsize));
    clear("params","F","yout","tout");
end

function params=defaultparams

    params.beta_N= 10^(1.7586);
    params.beta_G= 10^(1.8966);
    params.beta_F= 10^(0.0345);
    % params.beta_N= 1;
    % params.beta_G= 1;
    % params.beta_F= 1;

    params.mu = 1;
    params.nu = 1;
    params.xi = 1;
    % params.mu = 10^(1.68);
    % params.nu = 10^(-0.368);
    % params.xi = 10^(0.26);
    params.k1 = 1;

    params.n1 = 4;
    params.n2 = 4;

    params.sigma=0.5;     % noise amplitude in initial conditions
    params.P=12;        % number of cells per column
    params.Q=12;        % number of columns - MUST BE EVEN

function [dy] = li(t,y,params)
    beta_N = params.beta_N;
    beta_G = params.beta_G;
    beta_F = params.beta_F;

    mu = params.mu;
    nu = params.nu;
    xi = params.xi;
    k1 = params.k1;

    n1 = params.n1;
    n2 = params.n2;
    connectivity = params.connectivity;

    k = size(connectivity,2);
    synNotch_rtTA = y(1:k);
    rtTA = y(k+1:2*k);
    Gal4 = y(2*k+1:3*k);
    surEGFP = y(3*k+1:4*k);

    synNotch_rtTA_n = connectivity*synNotch_rtTA;
    surEGFP_n = connectivity*surEGFP;

    d_synNotch_rtTA = nu*(beta_N - k1*synNotch_rtTA.*surEGFP_n - synNotch_rtTA);
    d_rtTA = xi*(k1*synNotch_rtTA.*surEGFP_n - rtTA);
    d_Gal4 = beta_G * rtTA.^n1./(1+rtTA.^n1) - Gal4;
    d_surEGFP = mu*(beta_F./(1+Gal4.^n2) - k1*surEGFP.*synNotch_rtTA_n - surEGFP);

    dy = [d_synNotch_rtTA; d_rtTA; d_Gal4; d_surEGFP];
    

function M=getconnectivityM(P,Q)

    k=P*Q;          % number of cells
    M=zeros(k,k);   % connectivity matrix
    w=1/6;          % weight for interactions

    % calculating the connectivity matrix
    for s=1:k
        kneighbor=findneighborhex(s,P,Q); 
        for r=1:6
            M(s,kneighbor(r))=w;
        end
    end

function y0=getIC(params,k)

    U = rand(k,1) - 1/2;
    % u_sysNotch_rtTA=rand(k,1) - 1/2;
    % U_rtTA=rand(k,1) - 1/2;
    % U_Gal4=rand(k,1) - 1/2;
    % U_surEGFP=rand(k,1) - 1/2;
    epsilon=1e-5;   % multiplicative factor of Delta initial condition
    sigma=params.sigma;      % noise amplitude in initial conditions

    % sysNotch_rtTA0=params.beta_N*(1:2:2*k).'/k; % initial synNotch-rtTA levels
    % rtTA0=(1:2:2*k).'/k;          % initial rtTA levels
    % Gal40=params.beta_G*(1:2:2*k).'/k;          % initial Gal4 levels
    % surEGFP0=params.beta_F*(1:2:2*k).'/k; % initial surEGFP levels

    sysNotch_rtTA0 = params.beta_N*ones(k,1);
    rtTA0 = zeros(k,1);
    Gal40 = zeros(k,1);
    surEGFP0 = epsilon*params.beta_G*(1+sigma*U);

    % sysNotch_rtTA0=zeros(k,1); % initial synNotch-rtTA levels
    % rtTA0=zeros(k,1);          % initial rtTA levels
    % Gal40=zeros(k,1);          % initial Gal4 levels
    % surEGFP0=zeros(k,1); % initial surEGFP levels

    y0 = [sysNotch_rtTA0; rtTA0; Gal40; surEGFP0];

function plot2cells(tout,yout,k)

    figure
    clf
    for i=1:2
        subplot(1,2,i)
        plot(tout,yout(:,i),'-r','linewidth',2)   % plot synNotch-rtTA level
        hold on
        plot(tout,yout(:,k+i),'-b','linewidth',2) % plot rtTA levels
        plot(tout,yout(:,2*k+i),'-g','linewidth',2) % plot Gal4 levels
        plot(tout,yout(:,3*k+i),'-m','linewidth',2) % plot surEGFP levels
        title(['cell #',num2str(i)])
        xlabel('t [a.u]'); ylabel('concentration [a.u]')
        legend('synNotch-rtTA','rtTA','Gal4','surEGFP')
        hold off
    end

function out = findneighborhex(ind,P,Q)

    % This function finds the 6 neighbors of cell ind
    [p,q] = ind2pq(ind,P);

    % above and below:
    out(1) = pq2ind(mod(p,P)+1,q,P);
    out(2) = pq2ind(mod(p-2,P)+1,q,P);

    % left and right sides:
    qleft = mod(q-2,Q)+1;
    qright = mod(q,Q)+1;

    if q/2~=round(q/2),
        pup = p;
        pdown = mod(p-2,P)+1;
    else 
        pup = mod(p,P)+1;
        pdown = p;
    end;
    out(3) = pq2ind(pup,qleft,P);
    out(4) = pq2ind(pdown,qleft,P);
    out(5) = pq2ind(pup,qright,P);
    out(6) = pq2ind(pdown,qright,P);

function ind=pq2ind(p,q, P)
    ind = p + (q-1)*P;

function [p,q]=ind2pq(ind, P)
    q = 1+floor((ind-1)/P);
    p = ind - (q-1)*P;

function plotHexagon(p0,q0,c)

    % this function plots a hexagon centered at coordinates p,q

    s32 = sqrt(3)/4;
    q = q0*3/4;
    p = p0*2*s32;
    if q0/2 == round(q0/2),
    p = p+s32;
    end;

    x(1)=q-.5; x(2)=q-.25; x(3)=q+.25; 
    x(4)=q+.5; x(5)=q+.25; x(6)=q-.25;

    y(1)=p ; y(2)=p+s32; y(3)=p+s32; 
    y(4)=p; y(5)=p-s32; y(6)=p-s32;

    patch(x,y,c,'linewidth',2);

function F=movielattice(tout,yout,P,Q,k)

    % This function generates a movie of patterning in hexagonal 
    % lattice. The color represents the level of surEGFP. It also 
    % saves the movie as an AVI file. 
    figure
    clf
    Cmax=max(yout(end,1:k)); % finds max(surEGFP) at the end point
    Cmin=min(yout(end,1:k)); % finds min(surEGFP) at the end point
    frameind=0;
    F = moviein(numel(1:20:length(tout))); % preallocate movie structure
    for tind = length(tout):20:length(tout),   % shows every 5th frame
        clf;
        for i = 1:P,
            for j = 1:Q,
                ind = pq2ind(i,j,P);
                mycolor = min([yout(tind,ind)/Cmax,1]); 
                plotHexagon(i,j,[1-mycolor,1-mycolor,1]);
            end;
        end;
        axis image; axis off; box off;
        title(['Cmax = ',num2str(Cmax), '    ', 'Cmin = ',num2str(Cmin)]);
        frameind=frameind+1;
        F(frameind) = getframe; % generates a movie variable
    end;    
    % save movie in avi format
    video = VideoWriter('movielattice','MPEG-4');
    open(video);
    writeVideo(video, F);
    close(video);

function level = otsu(histogramCounts,n)
    total = sum(histogramCounts); % total number of pixels in the image 
    %% OTSU automatic thresholding
    top = n;
    sumB = 0;
    wB = 0;
    maximum = 0.0;
    sum1 = dot(0:top-1, histogramCounts);
    for ii = 1:top
        wF = total - wB;
        if wB > 0 && wF > 0
            mF = (sum1 - sumB) / wF;
            val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
            if ( val >= maximum )
                level = ii;
                maximum = val;
            end
        end
        wB = wB + histogramCounts(ii);
        sumB = sumB + (ii-1) * histogramCounts(ii);
    end
    if ~exist('level','var') || isnan(level)
        level = -1;
    end

function H=getCriterion(yout,P,Q,k,n)

    M=getconnectivityM(P,Q);

    yout_normlized=(yout-min(yout))/(max(yout)-min(yout)+1e-5);
    yout_integer = floor(yout_normlized*(n-1)+1);
    edges=0.5:1:(n+0.5);
    yout_hist=histcounts(yout_integer,edges);
    yout_threshold=otsu(yout_hist,n);

    H = zeros(numel(yout),2);
    for i=1:k
        H(i,2)=yout_integer(i)<yout_threshold;
        hex_n=find(M(i,:)~=0);
        for r=1:6
            H(i,1)=H(i,1)+double(yout_integer(hex_n(r))>=yout_threshold);
        end
    end
    H=H(H(:,2)~=0,1);
    H=median(H);


function yout_f = remove_padding(path_result, path_background)
    load(path_result);
    yout_f_result = yout_f;
    tout_result = tout;
    load(path_background);
    yout_f_background = yout_f;
    tout_background = tout;

    yout_f_result = normalize_yout(yout_f_result);
    yout_f_background = normalize_yout(yout_f_background);

    yout_f = yout_f_result - yout_f_background;
    yout_f = max(yout_f, 0);

function normalized = normalize_yout(yout)
    Cmax=max(yout(end,:)); % finds max(surEGFP) at the end point
    normalized = min(yout(end,:)/Cmax, 1);