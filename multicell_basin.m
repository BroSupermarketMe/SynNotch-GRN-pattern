% multicell_basin：统计不同细胞代数下，随机初始条件、随机几何条件（细胞生长分裂位置）的CV值
% i：细胞代数
for i = 18:2:20
    multicell_basin_i(i);
end

function multicell_basin_i(cell_generation)
    arguments
        cell_generation = 8;
    end
    simu_num = 1000;
    Hs_surEGFP = zeros(simu_num,6);
    Hs_Gal4 = zeros(simu_num,6);
    CVs_surEGFP = zeros(simu_num,1);
    CVs_Gal4 = zeros(simu_num,1);
    % 判断是否存在保存路径文件夹
    path_save = ['./figures_random_growth/cell_generation_' num2str(cell_generation)];
    if ~exist(path_save,'dir')
        mkdir(path_save);
    end
    parfor i = 1:simu_num
    % for i = 1:simu_num
        h_i = figure('position',[0,0,1800,900]);
        h_i.Visible = 'off';
        [~,~,~,cv_surEGFP,cv_Gal4,H_surEGFP,H_Gal4,isstable]=multicell_growth(cell_generation);
        exportgraphics(h_i,[path_save '/simulation_' num2str(i) '.png']);
        close(h_i);
        if isstable
            Hs_surEGFP(i,:) = H_surEGFP;
            Hs_Gal4(i,:) = H_Gal4;
            CVs_surEGFP(i) = cv_surEGFP;
            CVs_Gal4(i) = cv_Gal4;     
        end
    end
    save([path_save '/Hs_' num2str(simu_num) '.mat'],'Hs_surEGFP','Hs_Gal4','CVs_surEGFP','CVs_Gal4','-mat');
end