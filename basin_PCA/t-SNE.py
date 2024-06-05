import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans

def t_SNE(X):
    # 1. 初始化t-SNE模型
    tsne = TSNE(n_components=2, random_state=0)

    # 2. 适配t-SNE模型到数据上
    X_tsne = tsne.fit_transform(X)

    # 3. 返回降维后的特征矩阵
    return X_tsne

# 1. 读取.mat文件
mat_path = "C:\\Users\\JiajunTang\\Desktop\\20231106代码\\Lateral_inhibition-master\\figures\\cell_generation_10\\Hs.mat"
mat = scipy.io.loadmat(mat_path)

# 获取变量H_G和H_F
H_G = mat['Hs_Gal4']  # 从.mat文件中获取H_G
H_F = mat['Hs_surEGFP']  # 从.mat文件中获取H_F

# 过滤均为0的行
H_G = H_G[~np.all(H_G == 0, axis=1)]
H_F = H_F[~np.all(H_F == 0, axis=1)]

# 2. 转换为numpy数组
H_G = np.array(H_G)  # 转换为numpy数组
H_F = np.array(H_F)  # 转换为numpy数组

H = np.hstack((H_G, H_F))

# 3. t-SNE降维
H_tsne = t_SNE(H)

# 4. 绘制t-SNE降维后的散点图，H_G为红色，H_F为蓝色，两个subplots
plt.scatter(H_tsne[:, 0], H_tsne[:, 1], c='red', label='Gal4 & GFP concat')
plt.title('t-SNE Results of Gal4 & GFP concat')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend()

plt.show()