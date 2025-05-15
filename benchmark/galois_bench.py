import numpy as np
import galois
import json

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
    for _ in range(num_samples):
        rand_mat = GF.Random((n_rows, n_cols))
        row_reduced = rand_mat.row_reduce()
        input = rand_mat.flatten()
        output = row_reduced.flatten()
        rank = np.linalg.matrix_rank(row_reduced)
        data.append({
            "in": [int(input[ix]) for ix in range(len(input))],
            "out": [int(output[ix]) for ix in range(len(output))],
            "rank": rank,
        })
    json_dump["input"] = input_data
    json_dump["output"] = output_data
    with open("/Users/matt/repos/qec/benchmark/galois_bench_data.json", 'w') as f:
        json.dump({
            "num_samples": num_samples,
            "n_rows": n_rows,
            "n_cols": n_cols,
            "finite_field": q,
            "memory_layout": "RowMajor",
            "data": data
            }, f)