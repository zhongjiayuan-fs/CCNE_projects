# ccne_main.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# load entropy_matrix
mat = loadmat('pericyte_to_neuron.mat')
entropy = mat['entropy_matrix']    # shape: (genes, 623, 6)

cell_num   = [76, 86, 48, 283, 61, 69]
time_point = 6
top_ratio  = 0.05
total_nodes = entropy.shape[0]

SH = [ ]    
for t in range(time_point):
    te = entropy[:, :cell_num[t], t]           
    SH_t = [
        np.sum( np.sort(te[:,c])[::-1][:int(total_nodes*top_ratio)] )
        for c in range(cell_num[t])
    ]
    SH.append(SH_t)

SH = np.array([
    np.pad(SH[t], (0, max(cell_num)-len(SH[t])), 'constant')
    for t in range(time_point)
])


result = np.array([ np.mean(SH[t,:cell_num[t]]) for t in range(time_point) ])


data_by_group = [ SH[t, :cell_num[t]] for t in range(time_point) ]
positions     = list(range(2, 8))   
labels        = positions           


plt.figure(figsize=(8,6))


boxprops    = dict(color=(0.9,0.1,0.1))
medianprops = dict(color=(0.9,0.1,0.1))
plt.boxplot(
    data_by_group,
    notch=True, sym='r+',
    boxprops=boxprops,
    medianprops=medianprops,
    positions=positions
)


plt.plot(positions, result, color=(0.9,0.1,0.1), linewidth=3)
plt.scatter(positions, result, s=100, color=(0.9,0.1,0.1), zorder=3)


plt.ylabel('CSCNE')
plt.ylim(0, 40)
plt.yticks(np.arange(0,41,10))
plt.xticks(labels)
plt.box(False)
plt.tight_layout()
plt.show()
