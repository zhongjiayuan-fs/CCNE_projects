import numpy as np
from scipy.spatial.distance import cdist

def Rewiring_cnet(x, y, k=3):
    """
    DMI: compute the direct entropy causality from x to y.

    Parameters
    ----------
    x : array-like, shape (N,)
        N samples.
    y : array-like, shape (N,)
        N samples.
    k : int, optional (default=3)
        The k-th nearest neighbor to use in calculating entropy, at least 2.

    Returns
    -------
    XY_edge : ndarray, shape (N,)
        Direct entropy-causality value for each sample i (x_i → y_i).
    """
    # Ensure column-vector shape (N, 1)
    x = np.asarray(x).reshape(-1, 1)
    y = np.asarray(y).reshape(-1, 1)
    N = x.shape[0]

    # Pairwise Euclidean distances
    dx = cdist(x, x, metric='euclidean')
    dy = cdist(y, y, metric='euclidean')

    # k-NN indices (include self → k+1)
    idx = np.argsort(dx, axis=1)
    idy = np.argsort(dy, axis=1)

    XNNid = idx[:, :k + 1].copy()
    YNNid = idy[:, :k + 1].copy()

    # Force first column to be self indices (1:N in MATLAB → 0…N-1 here)
    XNNid[:, 0] = np.arange(N)
    YNNid[:, 0] = np.arange(N)

    # Gather y-values of x-NN 和 y-NN
    YxNN = y[XNNid]          # shape (N, k+1)
    YNN  = y[YNNid]          # shape (N, k+1)

    # ───────── YxNN Δ ─────────
    xfirst_column = YxNN[:, 0]  
    YxNN_dis = np.abs(YxNN[:, 1:] - xfirst_column[:, None])  
    YxNN_dis = YxNN_dis.squeeze()  

    # ───────── YNN Δ& epsilon ─────────
    yfirst_column = YNN[:, 0]
    YNN_dis = np.abs(YNN[:, 1:] - yfirst_column[:, None])
    YNN_dis = YNN_dis.squeeze()
    halfepsilon = np.max(YNN_dis, axis=1).squeeze()  # (N,)

    
    common_matrix = (YxNN_dis < halfepsilon[:, None])  # (N, k)
    common_num = np.sum(common_matrix, axis=1)  # (N,)

    c = common_num / N  # shape: (N,)
    a = (k * np.ones(N)) / N  # shape: (N,)
    b = (k * np.ones(N)) / N  # shape: (N,)
    XY_edge = c * np.log(c / (a * b))  # shape: (N,)

    XY_edge[XY_edge < 0] = 0
    XY_edge = np.nan_to_num(XY_edge)

    return XY_edge
