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

% paramsearch_sysNotch：参数扫描：在周期性边界条件下（multicell_sysNotch），
% 使用不同参数，获取稳态后的各项指标（Dmax/Dmin，达到稳态的时间，Gal4和GFP的H值
% （非活性细胞周围有多少个活性细胞）），并保存结果到result/paramsearch_sysNotch_result_new.mat

function paramsearch_sysNotch

    % This code plots the log(Dmax/Dmin) and the time required for patterning 
    % for different betaD and betaR. TO RUN FASTER, COMMENT OUT PLOT2CELLS AND 
    % MOVIELATTICE within multicell_LI function. For a better exploration in 
    % the parameter space, it is important to set a large enough Tmax parameter 
    % in multicell_LI function.
    
    
    % fixed parameters
    beta_N= logspace(-1,3,30);
    beta_G= logspace(-1,3,30);
    beta_F= logspace(-1,3,30);
    % params.beta_N= 10^(0.9655);
    % params.beta_G= 10^(1.2759);
    % params.beta_F= 10^(0.345);

    params.mu = 1;
    params.nu = 1;
    params.xi = 1;
    % mu = logspace(-1,2,20);
    % nu = logspace(-1,2,20);
    % xi = logspace(-1,2,20);
    params.k1 = 1;

    params.n1 = 4;
    params.n2 = 4;

    params.sigma=0.2;   % noise amplitude in initial conditions
    params.P=12;        % number of cells per column
    params.Q=12;        % number of columns - MUST BE EVEN

    k=params.Q*params.P;
    
    Dmax=zeros(length(beta_N),length(beta_G),length(beta_F));
    Dmin=zeros(length(beta_N),length(beta_G),length(beta_F));
    T=zeros(length(beta_N),length(beta_G),length(beta_F));
    H1=zeros(length(beta_N),length(beta_G),length(beta_F));
    H2=zeros(length(beta_N),length(beta_G),length(beta_F));
    % Dmax=zeros(length(mu),length(nu),length(xi));
    % Dmin=zeros(length(mu),length(nu),length(xi));
    % T=zeros(length(mu),length(nu),length(xi));
    % H1=zeros(length(mu),length(nu),length(xi));
    % H2=zeros(length(mu),length(nu),length(xi));

    iteration=(length(beta_N)*length(beta_G)*length(beta_F));
    % iteration=(length(mu)*length(nu)*length(xi));
    % h=ParforProgressbar(iteration); % generates a waitbar
    for ind=1:iteration

        [i,j,m]=ind2sub([length(beta_N),length(beta_G),length(beta_F)],ind);
        % [i,j,m]=ind2sub([length(mu),length(nu),length(xi)],ind);
        
        param = params;
        % param.mu = mu(i);
        % param.nu = nu(j);
        % param.xi = xi(m);
        param.beta_N = beta_N(i);
        param.beta_G = beta_G(j);
        param.beta_F = beta_F(m);

        % h.increment();
        [yout,tout,H0] = multicell_sysNotch(param); % calling the LI solver
        
        % finding max and min values of D
        Dmax(ind)=max(max(yout(end,2*k+1:3*k)));
        Dmin(ind)=abs(min(min(yout(end,2*k+1:3*k))));
        H1(ind)=H0(1);
        H2(ind)=H0(2);

        % finding cases where patterning occurs (when Dmax/Dmin>1.2)
        % and getting the patterning time 
        if Dmax(ind)/Dmin(ind)>1.2 
            T(ind)=getPatterningTime(tout,yout,k,Dmax(ind),Dmin(ind));
        else
            T(ind)=NaN; % patterning time is not set for the no patterning case
        end
    end
    % delete(h)
    % figure(23)
    % imagesc(log10(betaD),log10(betaR),log10(Dmax./Dmin));
    % set(gca,'YDir','normal')
    % xlabel('log(\beta_r)','fontsize',14);
    % ylabel('log(\beta_d)','fontsize',14);
    % title('log(d_{max}/d_{min})','fontsize',14)
    % colorbar
    
    % figure(24)
    % imagesc(log10(betaD),log10(betaR),T);
    % set(gca,'YDir','normal')
    % xlabel('log(\beta_r)','fontsize',14);
    % ylabel('log(\beta_d)','fontsize',14);
    % title('Time for patterning [normalized time]','fontsize',14)
    % colorbar
    save('result/paramsearch_sysNotch_result_new.mat','Dmax','Dmin','T','H1','H2','-mat')
    
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