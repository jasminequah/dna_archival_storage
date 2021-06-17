import numpy as np
import timeit
import torch
import matplotlib.pyplot as plt

def measure_sparse(hsize=3072, wsize=768, sparsities=None, gpu=False):
    print("size = (%s, %s)" % (hsize, wsize))
    if sparsities is None:
        sparsities = np.concatenate(np.arange(0., 0.9, 0.05), np.arange(0.9, 1., 0.005))
    times = np.zeros_like(sparsities, dtype=float)

    for i, sparsity in enumerate(sparsities):
        nnz = int((1.-sparsity) * hsize * wsize)
        # indices are coordinates of non-zero values in matrix
        setup = (
            'import torch;'
            'hsize = {0};'
            'wsize = {1};'
            'nnz = {2};'
            'device = torch.device("cuda:0" if {3} else "cpu");'
            'spmat = torch.sparse_coo_tensor(indices = torch.tensor([torch.randint(hsize,(nnz,)).tolist(),torch.randint(wsize,(nnz,)).tolist()]), values = torch.randn(nnz), dtype=torch.float32, device=device).coalesce();'
            'multmat = torch.randn(wsize,wsize).to(device)'
            ).format(hsize, wsize, nnz, gpu)
        reps, time = timeit.Timer(stmt = 'torch.sparse.mm(spmat,multmat)', setup = setup).autorange()
        times[i] = time/reps
        print("sparsity %.3f: %s" % (sparsity, time/reps))
    return times

def measure_dense(hsize=3072, wsize=768, sparsities=None, gpu=False):
    print("size = (%s, %s)" % (hsize, wsize))
    if sparsities is None:
        sparsities = np.concatenate(np.arange(0., 0.9, 0.05), np.arange(0.9, 1., 0.005))
    times = np.zeros_like(sparsities, dtype=float)

    for i,sparsity in enumerate(sparsities):
        nnz = int((1-sparsity) * hsize * wsize)
        # indices are coordinates of non-zero values in matrix
        setup = (
            'import torch;'
            'hsize = {0};'
            'wsize = {1};'
            'nnz = {2};'
            'device = torch.device("cuda:0" if {3} else "cpu");'
            'mat1 = torch.randn(hsize,wsize).to(device);'
            'mat2 = torch.randn(wsize,wsize).to(device)'
            ).format(hsize, wsize, nnz, gpu)
        reps, time = timeit.Timer(stmt = 'torch.mm(mat1,mat2)', setup = setup).autorange()
        times[i] = time/reps
        print("sparsity %.3f: %s" % (sparsity, time/reps))
    return times

sparse_sizes = [(3072, 768)]
sparsities = np.concatenate((np.linspace(0.75, 0.95, num=5, endpoint=False), np.linspace(0.95, 1., 30, endpoint=False)))
for sparse_size in sparse_sizes:
    hsize, wsize = sparse_size
    results = dict()
    results['sparse_cpu_time'] = measure_sparse(hsize, wsize, sparsities)
    results['dense_cpu_time'] = measure_dense(hsize, wsize, sparsities)
    if torch.cuda.is_available():
        results['sparse_gpu_time'] = measure_sparse(hsize, wsize, sparsities, gpu=True)
        results['dense_gpu_time'] = measure_dense(hsize, wsize, sparsities, gpu=True)

    plt.figure()
    plt.title("Matrix multiplication of (%sx%s) by (%sx%s)" % (hsize, wsize, wsize, wsize))
    for label, result in results.items():
        plt.plot(sparsities, result, label=label)
    plt.xlabel('Sparsity (%) of non-zero elements')
    plt.ylabel('Log-time (s)')
    plt.yscale("log")
    plt.legend()
    plt.savefig("sparse_vs_dense_%s_%s_cpu_log_new.png" % (hsize, wsize), dpi=800)
