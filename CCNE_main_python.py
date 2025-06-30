# ccne_main.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load entropy_matrix
mat = loadmat('pericyte_to_neuron.mat')
entropy = mat['entropy_matrix']    # shape: (genes, 623, 6)

cell_num   = [76, 86, 48, 283, 61, 69]
time_point = 6
top_ratio  = 0.05
total_nodes = entropy.shape[0]

# Calculate SH and result
SH = [ ]    
for t in range(time_point):
    te = entropy[:, :cell_num[t], t]           # genes × cells_at_t
    # For each cell 𝑐 genes are ranked and the top 5% are summed.
    SH_t = [
        np.sum( np.sort(te[:,c])[::-1][:int(total_nodes*top_ratio)] )
        for c in range(cell_num[t])
    ]
    SH.append(SH_t)

SH = np.array([
    np.pad(SH[t], (0, max(cell_num)-len(SH[t])), 'constant')
    for t in range(time_point)
])

# 每个时间点的平均 CSCNE
result = np.array([ np.mean(SH[t,:cell_num[t]]) for t in range(time_point) ])

# 组合数据 & 分组标签
data_by_group = [ SH[t, :cell_num[t]] for t in range(time_point) ]
positions     = list(range(2, 8))   # x=2,3,4,5,6,7
labels        = positions           # 直接用 2…7 作为 x 轴标签


plt.figure(figsize=(8,6))

# 箱线图
boxprops    = dict(color=(0.9,0.1,0.1))
medianprops = dict(color=(0.9,0.1,0.1))
plt.boxplot(
    data_by_group,
    notch=True, sym='r+',
    boxprops=boxprops,
    medianprops=medianprops,
    positions=positions
)

# 折线 + 散点
plt.plot(positions, result, color=(0.9,0.1,0.1), linewidth=3)
plt.scatter(positions, result, s=100, color=(0.9,0.1,0.1), zorder=3)

# 美化
plt.ylabel('CSCNE')
plt.ylim(0, 40)
plt.yticks(np.arange(0,41,10))
plt.xticks(labels)
plt.box(False)
plt.tight_layout()
plt.show()
