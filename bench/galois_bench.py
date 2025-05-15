import numpy as np
import galois
import json
import time

if __name__ == "__main__":
    num_samples = 200
    n_rows = 50
    n_cols = 70
    q = 5
    json_dump = {
        "num_samples": num_samples,
        "n_rows": n_rows,
        "n_cols": n_cols,
        "finite_field": q,
        "memory_layout": "RowMajor"
        }
    GF = galois.GF(q)
    input_data = []
    output_data = []
    data = []
    tot_time = 0.0
    for _ in range(num_samples):
        rand_mat = GF.Random((n_rows, n_cols))
        start = time.time()
        row_reduced = rand_mat.row_reduce()
        end = time.time()
        tot_time += end - start
        input = rand_mat.flatten()
        output = row_reduced.flatten()
        rank = np.linalg.matrix_rank(row_reduced)
        data.append({
            "input": [int(input[ix]) for ix in range(len(input))],
            "output": [int(output[ix]) for ix in range(len(output))],
            "rank": rank,
        })
    print("tot time: ", tot_time)
    print("time per matrix: ", tot_time / num_samples)
    json_dump["input"] = input_data
    json_dump["output"] = output_data
    with open("/Users/matt/repos/qec/bench/galois_bench_data.json", 'w') as f:
        json.dump({
            "num_samples": num_samples,
            "n_rows": n_rows,
            "n_cols": n_cols,
            "finite_field": q,
            "memory_layout": "RowMajor",
            "time_per_matrix": tot_time / num_samples,
            "data": data
            }, f)