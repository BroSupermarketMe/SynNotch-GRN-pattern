% multicell_pedigree：统计细胞自然生长下的CV值（细胞生长分裂位置随机）

% CVs_surEGFP = zeros(20,1);
% CVs_Gal4 = zeros(20,1);
% % load('./figures/CVs_pedigree.mat','CVs_surEGFP','CVs_Gal4');

% for i = 8:2:16
%     [~,~,~,cv_surEGFP,cv_Gal4]=multicell_growing(i);
%     CVs_surEGFP(i) = cv_surEGFP;
%     CVs_Gal4(i) = cv_Gal4;
% end

% save('./figures_random_growth/CVs_pedigree_slow.mat','CVs_surEGFP','CVs_Gal4','-mat');

simu_num = 100;

for i = 18:2:20
    CVs_surEGFP = zeros(simu_num,1);
    CVs_Gal4 = zeros(simu_num,1);
    
    parfor j = 1:simu_num
        [~,~,~,cv_surEGFP,cv_Gal4]=multicell_growing(i);
        CVs_surEGFP(j) = cv_surEGFP;
        CVs_Gal4(j) = cv_Gal4;
    end
    save_path = ['./figures_random_growth/cell_generation_' num2str(i)];
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end
    save([save_path 'CVs_pedigree_' num2str(simu_num) '.mat'],'CVs_surEGFP','CVs_Gal4','-mat');
end