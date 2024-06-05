function [Cell_N_generations,entropy, surround_surEGFP,surround_Gal4,cv_surEGFP,cv_Gal4,H_surEGFP,H_Gal4,isstable] = multicell_growing(cell_generation, para)
arguments
    cell_generation (1,1) double = 8
    para struct = defaultparams(cell_generation)
end

global CurrentTime Coordinate_X Coordinate_Y BirthTime DivisionTime CellAge LiveorDeath SpaceFull ...
    synNotch_rtTA rtTA Gal4 surEGFP tspan
CurrentTime = 1;    % 细胞当前所处细胞周期，如果细胞死亡，该值停留在细胞死亡的那一刻
Coordinate_X = 2;
Coordinate_Y = 3;
BirthTime = 4;
DivisionTime = 5;
CellAge = 6;
LiveorDeath = 7;
SpaceFull = 8;

synNotch_rtTA = 9;
rtTA = 10;
% rtTA_dox = 11;
% dox = 12;
Gal4 = 11;
surEGFP = 12;
CharacteristicsNumber = 12;  % 表征一个细胞的特征数目

COLOR_SET = 64; % 绘图的颜色数目

tspan = para.tspan;
% 细胞从出生到分裂所需要的时间，对应Gene_regulatory_function_2D更新Cellcyclelength-1次，因为细胞出生占一个Cellcyclelength
Cellcyclelength = para.Cellcyclelength;
max_Cellcyclelength = para.max_Cellcyclelength;
N = para.N; % 总代数N
TotalCellNumber = 2^(N+1)-1;  % 总细胞数目
last_Cellcyclelength = (N+1)*max_Cellcyclelength - N*Cellcyclelength;
% last_Cellcyclelength = max_Cellcyclelength;
% TotalTime = 2*(N+1);
Cell_Matrix = zeros(TotalCellNumber,CharacteristicsNumber);  % 用一个大的矩阵记录所有细胞

% F(N) = struct('cdata',[],'colormap',[]);
% frame = 1;

%% 定义初始细胞状态
time = 0;
% 选定初始细胞生长的位置，设定初始细胞的状态
cell_index = 1;
start_cell_position_x = 0;
start_cell_position_y = 0;
Cell_Matrix(cell_index,CurrentTime) = time;
Cell_Matrix(cell_index,Coordinate_X) = start_cell_position_x;
Cell_Matrix(cell_index,Coordinate_Y) = start_cell_position_y;
Cell_Matrix(cell_index,BirthTime) = time;
Cell_Matrix(cell_index,LiveorDeath) = 1;
Cell_Matrix(cell_index,CellAge) = time;
%  初始细胞内各基因表达量的初始状态
initial_synNotch_rtTA = 1;  % 参考那本书的设置
initial_rtTA = 0;
% initial_rtTA_dox = 0;
% initial_dox = 1;  % 目前dox设置为常数，一直维持初始值
initial_Gal4 = 0;
initial_surEGFP = 1e-5;
Cell_Matrix(cell_index,synNotch_rtTA) = initial_synNotch_rtTA;
Cell_Matrix(cell_index,rtTA) = initial_rtTA;
% Cell_Matrix(cell_index,rtTA_dox) = initial_rtTA_dox;
% Cell_Matrix(cell_index,dox) = initial_dox;
Cell_Matrix(cell_index,Gal4) = initial_Gal4;
Cell_Matrix(cell_index,surEGFP) = initial_surEGFP;

Cell_N_generations = cell(1,N+1);

%% 细胞周期
for time = 1:N
    %% 细胞生长
    [Cell_Matrix,i_generation_cell] = runcellcycle(Cell_Matrix, para, Cellcyclelength,false);
    Cell_N_generations(1,time) = Cell_Matrix(i_generation_cell);

    % 记录当前帧
    % [F,frame]=captureframe(Cell_Matrix,F,frame);

    %% 细胞分裂
    % 选择需要分裂的细胞，条件：1.存活；2.周围有空位
    divide_cell_number = find(Cell_Matrix(:,LiveorDeath)==1 & Cell_Matrix(:,SpaceFull) ~= 1);
    % s=RandStream('threefry4x64_20','Seed',123);
    s=RandStream('threefry4x64_20','Seed','shuffle');   % 随机数种子,随机生长
    randindex = randperm(s,numel(divide_cell_number));
    divide_cell_number_rand = divide_cell_number(randindex,:); %可以分裂的细胞分裂顺序随机
    for n_d = 1:size(divide_cell_number_rand,1)
        divide_cell_index = divide_cell_number_rand(n_d);
        daughter_cell_index_1 = 2*divide_cell_index;
        daughter_cell_index_2 = 2*divide_cell_index+1;
        neighborXYs = getneighborXY(Cell_Matrix(divide_cell_index,Coordinate_X:Coordinate_Y));
        
        % 检查细胞周围有没有空位（空位随分裂动态变化）
        % space_index = space_to_choose(Neighbor_Space(divide_cell_index,:));
        live_cell_inds=Cell_Matrix(:,LiveorDeath)==1;
        space_index = ~ismember(neighborXYs,Cell_Matrix(live_cell_inds,Coordinate_X:Coordinate_Y),'rows');
        if ~any(space_index)
            Cell_Matrix(divide_cell_index,SpaceFull)=1;
            continue
        end
        space_index = find(space_index);
        space_index = space_index(randperm(s,numel(space_index),1));
        neighbor_cell_coordinate_X = neighborXYs(space_index,1);
        neighbor_cell_coordinate_Y = neighborXYs(space_index,2);
        
        Cell_Matrix(daughter_cell_index_1,:) = Cell_Matrix(divide_cell_index,:);
        Cell_Matrix(daughter_cell_index_2,:) = Cell_Matrix(divide_cell_index,:);

        Cell_Matrix(daughter_cell_index_1,BirthTime) = time;
        Cell_Matrix(daughter_cell_index_1,CellAge) = 0;
        
        Cell_Matrix(daughter_cell_index_2,Coordinate_X) = neighbor_cell_coordinate_X;
        Cell_Matrix(daughter_cell_index_2,Coordinate_Y) = neighbor_cell_coordinate_Y;
        Cell_Matrix(daughter_cell_index_2,BirthTime) = time;
        Cell_Matrix(daughter_cell_index_2,CellAge) = 0;
            
        Cell_Matrix(divide_cell_index,DivisionTime) = time;
        Cell_Matrix(divide_cell_index,LiveorDeath) = 0;

        % 记录当前帧
        % [F,frame]=captureframe(Cell_Matrix,F,frame);
    end
end

s = RandStream('threefry4x64_20','Seed','shuffle');
%% 最后一个细胞周期
% cell_occupancy = find(Cell_Matrix(:,LiveorDeath)==1);
% Cell_Matrix(cell_occupancy,:) = getIC(para,Cell_Matrix(cell_occupancy,:),s);

% figure
[Cell_Matrix,cell_occupancy,~,connectivity,~,~,isstable] = ...
    runcellcycle(Cell_Matrix, para, last_Cellcyclelength,false,s);

Cell_N_generations(1,time+1) = Cell_Matrix(cell_occupancy);

% 获取next nearest邻居均不在边缘的细胞
cell_matrix_nn_core_index = get_nn_core(Cell_Matrix(cell_occupancy,SpaceFull),connectivity);

% 记录当前帧
% F=captureframe(Cell_Matrix,F,frame);

% 原来的选取细胞的方法
cell_compare = Cell_Matrix(:,LiveorDeath)==1 & Cell_Matrix(:,SpaceFull)==1;
cell_matrix_compare = Cell_Matrix(cell_compare,:);
% TODO
M_compare_index = cell_compare(cell_occupancy);
connectivity_compare = connectivity(M_compare_index,M_compare_index);

% 更新后的选取细胞的方法，不包括nn邻居为边缘细胞
cell_draw = cell_occupancy(cell_matrix_nn_core_index);
cell_matrix_draw = Cell_Matrix(cell_draw,:);
connectivity_draw = connectivity(cell_matrix_nn_core_index,cell_matrix_nn_core_index);


%% 计算灰度共生矩阵
[entropy,surround_surEGFP,surround_Gal4,cv_surEGFP,cv_Gal4,H_surEGFP, H_Gal4, threshold_surEGFP, threshold_Gal4,max_surEGFP,min_surEGFP,max_Gal4,min_Gal4]...
 = getcriterion(cell_matrix_draw,connectivity_draw,cell_matrix_compare,connectivity_compare,COLOR_SET);
% writerObj = VideoWriter('test.mp4','MPEG-4');
% writerObj.FrameRate=0.5;
% open(writerObj);
% for frame=1:size(F,2)
%     writeVideo(writerObj,F(frame));
% end
% close(writerObj);

% Cell_Matrix(cell_occupancy,Gal4) = norm_species(Cell_Matrix(cell_occupancy,Gal4),COLOR_SET,max_Gal4,min_Gal4);
% Cell_Matrix(cell_occupancy,surEGFP) = norm_species(Cell_Matrix(cell_occupancy,surEGFP),COLOR_SET,max_surEGFP,min_surEGFP);
% Cell_Matrix(cell_occupancy,Gal4) = Cell_Matrix(cell_occupancy,Gal4) >= threshold_Gal4;
% Cell_Matrix(cell_occupancy,surEGFP) = Cell_Matrix(cell_occupancy,surEGFP) >= threshold_surEGFP;

% figure
% t2 = tiledlayout(2,2);
% % subplot(2,2,1)
% nexttile(t2)
% plotlattice(Cell_Matrix(cell_occupancy,[Coordinate_X,Coordinate_Y,synNotch_rtTA]),COLOR_SET);
% title('synNotch-rtTA');
% % subplot(2,2,2)
% nexttile(t2)
% plotlattice(Cell_Matrix(cell_occupancy,[Coordinate_X,Coordinate_Y,rtTA]),COLOR_SET);
% title('rtTA');
% % subplot(2,2,3)
% nexttile(t2)
% plotlattice(Cell_Matrix(cell_occupancy,[Coordinate_X,Coordinate_Y,Gal4]),COLOR_SET);
% title('Gal4');
% % subplot(2,2,4)
% nexttile(t2)
% plotlattice(Cell_Matrix(cell_occupancy,[Coordinate_X,Coordinate_Y,surEGFP]),COLOR_SET);
% title('surEGFP');
% sgtitle('Final state')

% figure
% T2 = tiledlayout(1,2);
% nexttile(T2);
% H_edges = 1:7;
% histogram('BinEdges',H_edges,'BinCounts',H_surEGFP);
% title('surEGFP')
% nexttile(T2);
% histogram('BinEdges',H_edges,'BinCounts',H_Gal4);
% title('Gal4')

end

function [Cell_Matrix, cell_occupancy, neighbor_space,connectivity,yout,tout,isstable] = ...
    runcellcycle(Cell_Matrix, para, cycle_length, isplot, s)
    arguments
        Cell_Matrix
        para
        cycle_length
        isplot = false
        s = RandStream('threefry4x64_20','Seed','shuffle')
    end
    global CurrentTime Coordinate_X Coordinate_Y BirthTime DivisionTime CellAge LiveorDeath SpaceFull ...
        synNotch_rtTA rtTA Gal4 surEGFP tspan
    % 获取当前存活的细胞
    cell_occupancy = find(Cell_Matrix(:,LiveorDeath)==1);
    cell_matrix = Cell_Matrix(cell_occupancy,:);

    % 设置connectivity matrix k*k, neighbor_space k*6
    [connectivity,neighbor_space] = getconnectivityM(cell_matrix(:,Coordinate_X:Coordinate_Y));

    % 基因调控网络运行一个 cell cycle
    y0_synNotch_rtTA = cell_matrix(:,synNotch_rtTA);
    y0_rtTA = cell_matrix(:,rtTA);
    y0_Gal4 = cell_matrix(:,Gal4);
    y0_surEGFP = cell_matrix(:,surEGFP);
    y0 = vertcat(y0_synNotch_rtTA,y0_rtTA,y0_Gal4,y0_surEGFP);
    [tout,yout] = ode23(@Gene_regulatory_function_2D, cycle_length*tspan, y0,[],connectivity,para);
    
    k = size(cell_occupancy,1);
    cell_matrix(:,synNotch_rtTA) = yout(end,1:k);
    cell_matrix(:,rtTA) = yout(end,k+1:2*k);
    cell_matrix(:,Gal4) = yout(end,2*k+1:3*k);
    cell_matrix(:,surEGFP) = yout(end,3*k+1:4*k);
    % 更新细胞是否有空位
    cell_matrix(:,SpaceFull) = ~any(neighbor_space,2);
    % 更新细胞当前所处的周期和年龄
    cell_matrix(:,CurrentTime) = cell_matrix(:,CurrentTime)+1;
    cell_matrix(:,CellAge) = cell_matrix(:,CellAge)+1;

    Cell_Matrix(cell_occupancy,:) = cell_matrix;

    if isplot
        i=randi(s,k);
        plot(tout,yout(:,i),'-r','linewidth',2)   % plot synNotch-rtTA levels 
        hold on
        plot(tout,yout(:,k+i),'-b','linewidth',2) % plot rtTA levels
        plot(tout,yout(:,2*k+i),'-g','linewidth',2) % plot Gal4 levels
        plot(tout,yout(:,3*k+i),'-m','linewidth',2) % plot surEGFP levels
        title(['cell #',num2str(i)])
        xlabel('t [a.u]'); ylabel('concentration [a.u]')
        legend('synNotch-rtTA','rtTA','Gal4','surEGFP')
    end
    % 判断是否达到稳态
    T_surEGFP = getPatterningTime(tout,yout(:,3*k+1:4*k));
    T_Gal4 = getPatterningTime(tout,yout(:,2*k+1:3*k));
    T_max = para.N*para.tspan(2);
    if (T_surEGFP && T_Gal4) && (T_surEGFP < 0.8*T_max && T_Gal4 < 0.8*T_max)
        isstable = true;
    else
        isstable = false;
    end

end

function [dy] = Gene_regulatory_function_2D(t,y,connectivity,params)
    beta_N = params.beta_N;
    beta_G = params.beta_G;
    beta_F = params.beta_F;

    mu = params.mu;
    nu = params.nu;
    xi = params.xi;
    k1 = params.k1;

    n1 = params.n1;
    n2 = params.n2;

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
end

function [connectivity,neighbor_space] = getconnectivityM(coord_matrix)
    k = size(coord_matrix,1);
    connectivity = zeros(k,k);
    neighbor_space = false(k,6);
    
    for i=1:k
        neighbor_cell_xy = getneighborXY(coord_matrix(i,:));
        [neighbor_index, space_index] = ismember(coord_matrix,neighbor_cell_xy,'rows');
        connectivity(i,:) = neighbor_index.';
        space_to_check = 1:6;
        space_index = ismember(space_to_check,space_index);
        neighbor_space(i,:) = ~space_index;
    end
end

function [neighborXY] = getneighborXY(xy_k)
    x_k = xy_k(1);
    y_k = xy_k(2);

    vertup = mod(x_k,2);
    vertdown = mod(x_k +1,2);
    neighborXY = [x_k,  y_k+1
                  x_k+1,y_k+vertup
                  x_k+1,y_k-vertdown
                  x_k,  y_k-1
                  x_k-1,y_k-vertdown
                  x_k-1,y_k+vertup];
end

function cell_matrix_nn_core_index = get_nn_core(space_matrix,connectivity)
    % space_matrix (k,1)
    k = size(space_matrix,1);
    cell_matrix_nn_core = zeros(k,6);
    for i=1:k
        if ~space_matrix(i)
            continue
        end
        n_indexs = find(connectivity(i,:));
        for j=1:6
            n_cell_index = n_indexs(j);
            cell_matrix_nn_core(i,j) = double(space_matrix(n_cell_index));
        end
    end
    cell_matrix_nn_core_index = sum(cell_matrix_nn_core,2)==6;
end

function [entropy,hlh_surEGFP,hlh_Gal4,cv_surEGFP,cv_Gal4,H_surEGFP,H_Gal4,threshold_surEGFP,threshold_Gal4,max_surEGFP,min_surEGFP,max_Gal4,min_Gal4]...
     = getcriterion(cell_matrix_draw, connectivity_draw,cell_matrix_compare,connectivity_compare,COLOR_SET)
    global CurrentTime Coordinate_X Coordinate_Y BirthTime DivisionTime CellAge LiveorDeath SpaceFull ...
        synNotch_rtTA rtTA Gal4 surEGFP
    % 三个角度：30,90,150度
    glcm_surEGFP = zeros(3,COLOR_SET,COLOR_SET);
    glcm_Gal4 = zeros(3,COLOR_SET,COLOR_SET);
    
    [cell_matrix_draw(:,surEGFP),max_surEGFP,min_surEGFP] = norm_species(cell_matrix_draw(:,surEGFP),COLOR_SET);
    [cell_matrix_draw(:,Gal4),max_Gal4,min_Gal4] = norm_species(cell_matrix_draw(:,Gal4),COLOR_SET);

    cell_matrix_compare(:,surEGFP) = norm_species(cell_matrix_compare(:,surEGFP),COLOR_SET,max_surEGFP,min_surEGFP);
    cell_matrix_compare(:,Gal4) = norm_species(cell_matrix_compare(:,Gal4),COLOR_SET,max_Gal4,min_Gal4);

    edges = 0.5:1:(COLOR_SET+0.5);
    histogram_surEGFP = histcounts(cell_matrix_draw(:,surEGFP),edges);
    threshold_surEGFP = otsu(histogram_surEGFP,COLOR_SET);
    histogram_Gal4 = histcounts(cell_matrix_draw(:,Gal4),edges);
    threshold_Gal4 = otsu(histogram_Gal4,COLOR_SET);
    surround_surEGFP = zeros(size(cell_matrix_draw,1),3);
    surround_Gal4 = zeros(size(cell_matrix_draw,1),3);

    for i = 1:size(cell_matrix_draw,1)
        cell_draw_coord = cell_matrix_draw(i,Coordinate_X:Coordinate_Y);
        cell_draw_surEGFP = cell_matrix_draw(i,surEGFP);
        cell_draw_Gal4 = cell_matrix_draw(i,Gal4);
        i_compare = cell_matrix_compare(:,Coordinate_X)==cell_draw_coord(1) & ...
            cell_matrix_compare(:,Coordinate_Y)==cell_draw_coord(2);
        neighbor_index = find(connectivity_compare(i_compare,:));
        
        % threshold==-1时，物质浓度是单峰，代表没有pattern
        surround_surEGFP(i,2) = (cell_draw_surEGFP<threshold_surEGFP)&(threshold_surEGFP>=0);
        surround_Gal4(i,2) = (cell_draw_Gal4<threshold_Gal4)&(threshold_Gal4>=0);

        surround_surEGFP(i,3) = (cell_draw_surEGFP>=threshold_surEGFP)&(threshold_surEGFP>=0);
        surround_Gal4(i,3) = (cell_draw_Gal4>=threshold_Gal4)&(threshold_Gal4>=0);

        for j=1:size(neighbor_index,2)
            neighbor_coord = cell_matrix_compare(neighbor_index(j),Coordinate_X:Coordinate_Y);
            neighbor_surEGFP = cell_matrix_compare(neighbor_index(j),surEGFP);
            neighbor_Gal4 = cell_matrix_compare(neighbor_index(j),Gal4);
            theta=getangle(cell_draw_coord, neighbor_coord);
            glcm_surEGFP(theta,cell_draw_surEGFP,neighbor_surEGFP) = ...
                glcm_surEGFP(theta,cell_draw_surEGFP,neighbor_surEGFP)+1;
            glcm_Gal4(theta,cell_draw_Gal4,neighbor_Gal4) = ...
                glcm_Gal4(theta,cell_draw_Gal4,neighbor_Gal4)+1;

            surround_surEGFP(i,1) = ...
                    surround_surEGFP(i,1) + double(neighbor_surEGFP>=threshold_surEGFP);
            surround_Gal4(i,1) = ...
                surround_Gal4(i,1) + double(neighbor_Gal4>=threshold_Gal4);   
        end
    end
    [entropy_surGFP, entropy_Gal4] = glcm_entropy(glcm_surEGFP, glcm_Gal4);
    entropy = harmmean([entropy_Gal4,entropy_surGFP]);

    H_edges = 1:7;
    hlh_surEGFP = surround_surEGFP(surround_surEGFP(:,2)~=0,1);
    hhh_surEGFP = surround_surEGFP(surround_surEGFP(:,3)~=0,1);
    H_surEGFP = histcounts(hlh_surEGFP,H_edges);
    hlh_Gal4 = surround_Gal4(surround_Gal4(:,2)~=0,1);
    hhh_Gal4 = surround_Gal4(surround_Gal4(:,3)~=0,1);
    H_Gal4 = histcounts(hlh_Gal4,H_edges);
    hlh_surEGFP = median(hlh_surEGFP,"all");
    hlh_Gal4 = median(hlh_Gal4,"all");

    cv_surEGFP = std(hhh_surEGFP)/mean(hhh_surEGFP);
    cv_Gal4 = std(hhh_Gal4)/mean(hhh_Gal4);

    % H_surGFP, H_Gal4需要正则化为概率密度
    H_surEGFP = H_surEGFP/sum(H_surEGFP);
    H_Gal4 = H_Gal4/sum(H_Gal4);

end

function theta = getangle(my_coord, neighbor_coord)
    neighborXY = getneighborXY(my_coord);
    [~,theta] = ismember(neighbor_coord,neighborXY,'rows');
    theta = mod(theta,3)+1;
end


function [entropy_surGFP, entropy_Gal4] = glcm_entropy(glcm_surEGFP, glcm_Gal4)
    % 输入为没有归一化的glcm，需要先分别归一化
    for i=1:size(glcm_surEGFP,1)
        glcm_surEGFP(i,:,:) = glcm_surEGFP(i,:,:) / sum(glcm_surEGFP(i,:,:), 'all');
        glcm_Gal4(i,:,:) = glcm_Gal4(i,:,:) / sum(glcm_Gal4(i,:,:), 'all');
    end
    entropy_surGFP = zeros(1,3);
    entropy_Gal4 = zeros(1,3);
    for i=1:size(glcm_surEGFP,1)
        entropy_surGFP(i) = -sum(glcm_surEGFP(i,:,:).*log10(glcm_surEGFP(i,:,:)),"all","omitnan");
        entropy_Gal4(i) = -sum(glcm_Gal4(i,:,:).*log10(glcm_Gal4(i,:,:)),"all","omitnan");
    end
end

function level = otsu(histogramCounts,COLOR_SET)
    total = sum(histogramCounts); % total number of pixels in the image 
    %% OTSU automatic thresholding
    top = COLOR_SET;
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
end

function [all_species,max_species,min_species] = norm_species(all_species,COLOR_SET,max_species,min_species)
    arguments
        all_species
        COLOR_SET
        max_species = max(all_species,[],'all');
        min_species = min(all_species,[],'all');
    end
    norm = max_species-min_species;
    if norm==0
        norm=0.00001;
    end 
    % 保证区间左开右闭
    all_species = max(min((all_species-min_species)/norm,1-1e-5),0);
    all_species = floor(all_species*COLOR_SET+1);
end


function plotHexagon(x0,y0,c)
    
    % this function plots a hexagon centered at coordinates p,q
    
    s32 = sqrt(3)/4;
    q = x0*3/4;
    p = y0*2*s32;
    if x0/2 ~= round(x0/2)
       p = p+s32;
    end
    
    x(1)=q-.5; x(2)=q-.25; x(3)=q+.25; 
    x(4)=q+.5; x(5)=q+.25; x(6)=q-.25;
    
    y(1)=p ; y(2)=p+s32; y(3)=p+s32; 
    y(4)=p; y(5)=p-s32; y(6)=p-s32;
    
    patch(x,y,c,'linewidth',2);
end

function plotlattice(cell_matrix,COLOR_SET)

% This function generates a movie of patterning in hexagonal 
% lattice. The color represents the level of Delta. It also 
% saves the movie as an AVI file.

% cell_matrix: (cell_index,3); dim=2: x,y,Gal4/GFP

cell_matrix(:,3) = norm_species(cell_matrix(:,3),COLOR_SET);

for i = 1:size(cell_matrix,1)
    x0=cell_matrix(i,1);
    y0=cell_matrix(i,2);
    mycolor=min([cell_matrix(i,3)/COLOR_SET,1]);
    plotHexagon(x0,y0,[1-mycolor,1-mycolor,1]);
end
axis image; axis off; box off; 
end

function [F,frame]=captureframe(Cell_Matrix,F,frame)
    global CurrentTime Coordinate_X Coordinate_Y BirthTime DivisionTime CellAge LiveorDeath SpaceFull ...
        synNotch_rtTA rtTA Gal4 surEGFP
    cell_occupancy = Cell_Matrix(:,LiveorDeath) ==1;
    plotlattice(Cell_Matrix(cell_occupancy,[Coordinate_X,Coordinate_Y,surEGFP]),9)
    title(sprintf('Cycle %d Division %d',time,n_d))
    F(frame) = getframe(gcf);
    frame = frame+1;
end


function params=defaultparams(cell_generation)
    params.tspan=[0 5];
    % params.Cellcyclelength = 0.4;
    params.Cellcyclelength = 2;
    params.max_Cellcyclelength = 4;
    params.N = cell_generation; % N>0

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

    params.sigma=1;     % noise amplitude in initial conditions
    params.P=12;        % number of cells per column
    params.Q=12;        % number of columns - MUST BE EVEN

end

function cell_matrix=getIC(params,cell_matrix,s)
    global synNotch_rtTA rtTA Gal4 surEGFP
    k = size(cell_matrix,1);

    % U = rand(k,1) - 1/2;
    U_sysNotch_rtTA=rand(s,k,1) - 1/2;
    U_rtTA=rand(s,k,1) - 1/2;
    U_Gal4=rand(s,k,1) - 1/2;
    U_surEGFP=rand(s,k,1) - 1/2;
    epsilon=1e-5;   % multiplicative factor of Delta initial condition
    sigma=params.sigma;      % noise amplitude in initial conditions

    % sysNotch_rtTA0=params.beta_N*(1:2:2*k).'/k; % initial synNotch-rtTA levels
    % rtTA0=(1:2:2*k).'/k;          % initial rtTA levels
    % Gal40=params.beta_G*(1:2:2*k).'/k;          % initial Gal4 levels
    % surEGFP0=params.beta_F*(1:2:2*k).'/k; % initial surEGFP levels

    sysNotch_rtTA0 = params.beta_N*(1+sigma * U_sysNotch_rtTA);
    rtTA0 = epsilon*(1+sigma*U_rtTA);
    Gal40 =epsilon*params.beta_G*(1+sigma*U_Gal4);
    surEGFP0 = epsilon*params.beta_F*(1+sigma*U_surEGFP);

    % sysNotch_rtTA0=zeros(k,1); % initial synNotch-rtTA levels
    % rtTA0=zeros(k,1);          % initial rtTA levels
    % Gal40=zeros(k,1);          % initial Gal4 levels
    % surEGFP0=zeros(k,1); % initial surEGFP levels

    cell_matrix(:,synNotch_rtTA) = sysNotch_rtTA0;
    cell_matrix(:,rtTA) = rtTA0;
    cell_matrix(:,Gal4) = Gal40;
    cell_matrix(:,surEGFP) = surEGFP0;
end

function T=getPatterningTime(tout,yout)
    
    % This function estimates patterning time. This is done by the follwoing 3 
    % steps:
    % 1. find all the high D cells ('onCells')
    % 2. find the time it takes for each 'on cell' to reach 90% of its final
    % level ('TonCells')
    % 3. get median value of the times calculated in stage 2

    % 在paramsearch_sysNotch基础上修改

    k = size(yout,2);
    Dmax = max(max(yout(end,:)));
    Dmin = abs(min(min(yout(end,:))));
    
    if Dmax/Dmin>1.2
        onCells=find(yout(end,:)>0.5*(Dmax+Dmin));
        for i=1:length(onCells)
            Tind=find(yout(:,onCells(i))>0.9*yout(end,onCells(i)),1,'first');
            TonCells(i)=tout(Tind);
        end 
        
        T=median(TonCells);
    else
        T=false;
    end
end