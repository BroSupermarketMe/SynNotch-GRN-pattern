import scipy.io
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt

from yellowbrick.cluster.elbow import kelbow_visualizer
from yellowbrick.datasets.loaders import load_nfl


def get_optimal_dim(data, name):
    # 初始化PCA模型
    pca = PCA()

    # 适配PCA模型到数据上
    pca.fit(data)  # 适配到H_G上，你也可以选择stacked_matrix

    # 获取最佳降维维度m
    explained_variance_ratio = pca.explained_variance_ratio_
    cumulative_explained_variance_ratio = np.cumsum(explained_variance_ratio)
    # 选择合适的维度m，使得累积解释方差比例达到一定的阈值，或者根据其他准则选择
    m = np.argmax(cumulative_explained_variance_ratio > 0.95) + 1  # 选择解释方差比例超过95%的维度

    # 使用选择的维度m来执行降维
    pca = PCA(n_components=m)
    data_pca = pca.fit_transform(data)
    # 如果选择了合并矩阵，则这里也需要相应地对合并后的矩阵执行降维

    # 现在，H_G_pca就是降维后的特征矩阵，形状为(n, m)，可以用于聚类分析
    print(f'shape of {name}_pac: ', data_pca.shape)
    return data_pca


def draw_cluster_plot(data_pca, name, N, n_clst=3):
    kmeans = KMeans(n_clusters=n_clst, random_state=0)
    labels = kmeans.fit_predict(data_pca)

    # 3. 可视化聚类结果（在前两个主成分上绘制散点图）
    plt.figure(figsize=(10, 8))

    # 计算每个数据点的密度估计
    kde = KernelDensity(bandwidth=0.5, metric='euclidean', kernel='gaussian')
    kde.fit(data_pca)
    log_density = kde.score_samples(data_pca)  # 计算对数密度估计

    # 根据密度值调整点的大小
    min_size = 10
    max_size = 100
    sizes = (log_density - np.min(log_density)) / (np.max(log_density) - np.min(log_density))  # 将密度值映射到 [0, 1] 范围
    sizes = min_size + sizes * (max_size - min_size)  # 映射到指定范围 [min_size, max_size]

    # 根据聚类结果绘制不同类别的数据点
    for label in range(n_clst):  # 根据聚类数目循环绘制
        plt.scatter(data_pca[labels == label, 0],  # x轴为第一个主成分
                    data_pca[labels == label, 1],  # y轴为第二个主成分
                    label=f'Cluster {label + 1}',  # 类别标签
                    alpha=0.8,  # 设置透明度
                    s=sizes[labels == label])  # 设置点的大小（根据密度估计值调整）

    # 绘制聚类中心
    centers = kmeans.cluster_centers_
    plt.scatter(centers[:, 0], centers[:, 1],  # 聚类中心的坐标
                c='red', s=200, alpha=1, label='Cluster Centers')  # 设置为红色圆圈表示

    # 设置图形属性
    plt.title(f'K-means Clustering Results of {name} (PCA: 2D)')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend()
    plt.grid(True)

    # 显示图形
    # plt.show()

    # 保存图片
    plt.savefig(f"C:\\Users\\JiajunTang\\Desktop\\20231106代码\\Lateral_inhibition-master\\figures\\cell_generation_{N}\\{name}.png")

    return labels


def find_outlier(data, labels, out_label):
    outlier_indexs = np.where(labels == out_label - 1)[0]
    # 随机选取indexs中的一个索引
    outlier_index = np.random.choice(outlier_indexs)
    outlier = data[outlier_index]
    return outlier


def main(N, CV_G_i, CV_F_i, CV_G_slow_i, CV_F_slow_i):
    # 1. 读取.mat文件
    mat_path = f"C:\\Users\\JiajunTang\\Desktop\\20231106代码\\Lateral_inhibition-master\\figures\\cell_generation_{N}\\Hs_200.mat"
    mat = scipy.io.loadmat(mat_path)

    # 获取变量H_G和H_F
    H_G = mat['Hs_Gal4']  # 从.mat文件中获取H_G
    H_F = mat['Hs_surEGFP']  # 从.mat文件中获取H_F
    # 变异系数
    CV_G = mat['CVs_Gal4']
    CV_F = mat['CVs_surEGFP']

    # 过滤均为0的行
    H_G = H_G[~np.all(H_G == 0, axis=1)]
    H_F = H_F[~np.all(H_F == 0, axis=1)]

    # 过滤CV_G, CV_F均为0的行，不过滤含NaN的行
    CV_G = CV_G[~np.all(CV_G == 0, axis=1)]
    CV_F = CV_F[~np.all(CV_F == 0, axis=1)]

    # 2. 转换为numpy数组
    H_G = np.array(H_G)  # 转换为numpy数组
    H_F = np.array(H_F)  # 转换为numpy数组

    # 3. 合并两个矩阵（如果需要）
    # 在这里你可以选择将H_G和H_F合并成一个大的矩阵，如果它们表示相似的数据类型
    stacked_matrix = np.hstack((H_G, H_F))  # 水平堆叠

    # 4. 执行PCA降维
    # H_G_pca = get_optimal_dim(H_G, 'H_G')
    # H_F_pca = get_optimal_dim(H_F, 'H_F')
    # H_pca = get_optimal_dim(stacked_matrix, 'H_G and H_F')
    pca = PCA(n_components=2)
    H_pca = pca.fit_transform(stacked_matrix)
    H_F_pca = pca.fit_transform(H_F)
    H_G_pca = pca.fit_transform(H_G)

    # 不执行降维
    # H_G_pca = H_G
    # H_F_pca = H_F

    # # 计算不同聚类个数的K-means模型的簇内平方和
    # inertias = []
    # for k in range(1, 10):  # 尝试不同的聚类个数
    #     kmeans = KMeans(n_clusters=k, random_state=0)
    #     kmeans.fit(H_G_pca)
    #     inertias.append(kmeans.inertia_)
    #
    # # 绘制肘部法则图像
    # plt.plot(range(1, 10), inertias, marker='o')
    # plt.xlabel('Number of Clusters')
    # plt.ylabel('Inertia')
    # plt.title('Elbow Method for Optimal K')
    # plt.show()

    # Use the quick method and immediately show the figure
    # kelbow_visualizer(KMeans(random_state=0), H_G_pca, k=(1,10))
    # kelbow_visualizer(KMeans(random_state=0), H_F_pca, k=(1,10))
    # kelbow_visualizer(KMeans(random_state=0), H_pca, k=(1,10))
    # 最佳都是4类
    # labels_G = draw_cluster_plot(H_G_pca, 'Gal4', N=N, n_clst=2)
    # labels_F = draw_cluster_plot(H_F_pca, 'surEGFP', N=N, n_clst=2)
    # labels_H = draw_cluster_plot(H_pca, 'Gal4 and surEGFP', N=N, n_clst=3)

    # outlier_G = find_outlier(H_G,labels_G,4)
    # outlier_F = find_outlier(H_F,labels_F,4)

    # 绘制outlier的直方图
    # plt.figure(figsize=(10, 8))
    # plt.hist(outlier_G, bins=6, label='Gal4', color='blue')
    # plt.show()
    # tangjiajun is pig

    plt.figure(figsize=(12, 8))
    plt.hist(CV_G, label='Gal4', color='blue', alpha=0.5)
    # 用对应颜色标注CV_G_i的位置
    plt.axvline(CV_G_i, color='blue', linewidth=2)
    plt.axvline(CV_G_slow_i, color='blue', linewidth=2, linestyle='--')
    plt.hist(CV_F, label='surEGFP', color='green', alpha=0.5)
    # 用对应颜色标注CV_F_i的位置
    plt.axvline(CV_F_i, color='green', linewidth=2)
    plt.axvline(CV_F_slow_i, color='green', linewidth=2, linestyle='--')
    plt.xlabel('CV')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(f"C:\\Users\\JiajunTang\\Desktop\\20231106代码\\Lateral_inhibition-master\\figures\\cell_generation_{N}\\CV.png")

CV_path = f"C:\\Users\\JiajunTang\\Desktop\\20231106代码\\Lateral_inhibition-master\\figures\\CVs_pedigree.mat"
CV_mat = scipy.io.loadmat(CV_path)
CV_slow_path = f"C:\\Users\\JiajunTang\\Desktop\\20231106代码\\Lateral_inhibition-master\\figures\\CVs_pedigree_slow.mat"
CV_slow_mat = scipy.io.loadmat(CV_slow_path)
CV_G = CV_mat['CVs_Gal4']
CV_F = CV_mat['CVs_surEGFP']
CV_G_slow = CV_slow_mat['CVs_Gal4']
CV_F_slow = CV_slow_mat['CVs_surEGFP']
for i in range(8, 21):
    CV_G_i = CV_G[i-1]
    CV_F_i = CV_F[i-1]
    CV_G_slow_i = CV_G_slow[i-1]
    CV_F_slow_i = CV_F_slow[i-1]
    main(i, CV_G_i, CV_F_i, CV_G_slow_i, CV_F_slow_i)