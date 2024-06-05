# SynNotch-GRN-pattern
Exploring the Impact of Cell Growth and Division on Pattern Formation: A Study Based on SynNotch Gene Regulatory Network
本科毕业设计《以SynNotch基因调控回路为例探究细胞生长分裂对模式形成的影响》代码

- multicell_basin：统计不同细胞代数下，随机初始条件、随机几何条件（细胞生长分裂位置）的CV值

- multicell_growth：用于模拟一个指定细胞代数N下的多细胞系统，不同初始条件后的稳态斑图，被multicell_basin循环使用

- multicell_pedigree：统计细胞自然生长下的CV值（细胞生长分裂位置随机）

- multicell_growing: 用于模拟一个自然生长的多细胞系统，被multicell_pedigree循环使用，与multicell_growth的不同在于每次分裂运行基因调控网络

- paramsearch_sysNotch：参数扫描：在周期性边界条件下（multicell_sysNotch），使用不同参数，获取稳态后的各项指标（Dmax/Dmin，达到稳态的时间，Gal4和GFP的H值（非活性细胞周围有多少个活性细胞）），并保存结果到result/paramsearch_sysNotch_result_new.mat

- multicell_sysNotch：模拟周期性边界条件下的多细胞斑图，被paramsearch_sysNotch循环使用

- screen_params: 参数扫描：在自然生长条件下（check_params），使用不同参数，获取稳态后的各项指标（Dmax/Dmin，达到稳态的时间，Gal4和GFP的H值（非活性细胞周围有多少个活性细胞）），并保存结果到MatData/result.mat以及MatData/params.mat

- check_params：模拟一个自然生长的多细胞系统，被screen_params循环使用

multicell_growth，multicell_growing，check_params主要代码基本上都差不多，check_params多了一些可视化的代码

- HFFT文件夹下的代码是傅里叶变换相关的

- basin_PCA文件夹下代码用于多细胞生长蒙特卡洛模拟的结果可视化
