import time
import numpy as np
import pandas as pd
from scipy.io import savemat
from Rewiring_cnet import Rewiring_cnet   
start_time = time.time()   


raw_df = pd.read_excel(
    'pericyte_to_neuron_data.xlsx',
    header=None,
    usecols=range(1, 624),
    engine='openpyxl'
)
raw_df = raw_df.apply(pd.to_numeric, errors='coerce')
mydata = raw_df.fillna(0.0).to_numpy(dtype=np.float64)
mydata = np.log(mydata + 1.0)  


adjacent_network = []
with open('Undirected_ppi_network.txt', 'r') as fid:
    for line in fid:
        parts = line.rstrip('\n').split('\t')
        adjacent_network.append(parts)
total_node_num = len(adjacent_network)

cell_num= [76, 86, 48, 283, 61, 69]   
time_point = 6
eps= 1e-10

cumsum_cell = np.cumsum(cell_num)
time_slices = [(0, cumsum_cell[0])] + [
    (cumsum_cell[i-1], cumsum_cell[i]) for i in range(1, time_point)
]


max_cells = max(cell_num)
entropy_matrix = np.zeros((total_node_num, max_cells, time_point), dtype=np.float64)


from joblib import Parallel, delayed

def process_center_gene(na, input_data, neighbor_list, k, N_cols, center_idx):
    weights = []
    for nei in neighbor_list:
        nei_idx = int(nei) - 1
        x = input_data[center_idx, :]
        y = input_data[nei_idx, :]
        out_x2y = Rewiring_cnet(x, y, k)
        weights.append(out_x2y)
    if not weights:
        return np.zeros(N_cols)
    local_w = np.vstack(weights)
    col_sums = local_w.sum(axis=0, keepdims=True)
    normalized = np.divide(
        local_w, col_sums, out=np.zeros_like(local_w), where=col_sums > 0)
    cols_ent = np.zeros(normalized.shape[1], dtype=np.float64)
    for col in range(normalized.shape[1]):
        p = normalized[:, col]
        nonzero = np.count_nonzero(p)
        center_data = input_data[center_idx, :].copy()
        sd1 = np.std(center_data, ddof=1)
        center_data = np.delete(center_data, col)
        sd2 = np.std(center_data, ddof=1)
        sd_delta = np.sqrt(N_cols) * abs(sd1 - sd2)
        H = -np.sum(p * np.log2(p + 1e-10))
        if nonzero > 0:
            cols_ent[col] = (1.0 / nonzero) * H * sd_delta
    return cols_ent

for t in range(time_point):
    start_col, end_col = time_slices[t]
    input_data = mydata[:, start_col:end_col]
    N_cols = input_data.shape[1]
    k = int(round(0.1 * N_cols))
    results = Parallel(n_jobs=12)(delayed(process_center_gene)(
        na, input_data, adjacent_network[na][1:], k, N_cols, int(adjacent_network[na][0])-1
    ) for na in range(total_node_num))
    for na in range(total_node_num):
        entropy_matrix[na, :N_cols, t] = results[na]
    print(f"Finished time point {t+1}/{time_point}")


entropy_matrix = np.nan_to_num(
    entropy_matrix,
    nan=0.0,
    posinf=0.0,
    neginf=0.0
)
entropy_matrix[entropy_matrix < 0] = 0.0

savemat('pericyte_to_neuron.mat', {'entropy_matrix': entropy_matrix})
end_time = time.time()   
print(f"Total runtime of the code：{end_time - start_time:.2f} s")