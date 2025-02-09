import numpy as np
import galois
import json

def create_benchmark_matrices():
    num_samples = 200
    n_rows = 50
    n_cols = 75
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
    for _ in range(num_samples):
        rand_mat = GF.Random((n_rows, n_cols), seed=1)
        row_reduced = rand_mat.row_reduce()
        input = rand_mat.reshape((n_rows * n_cols, 1))
        output = row_reduced.reshape((n_rows * n_cols, 1))
        for ix in range(len(input)):
            input_data.append(int(input[ix]))
        for ix in range(len(output)):
            output_data.append(int(output[ix]))
    json_dump["input"] = input_data
    json_dump["output"] = output_data
    with open("/Users/matt/repos/qec/scripts/galois_bench_data.json", 'w') as f:
        json.dump(json_dump, f)

if __name__ == "__main__":
    create_benchmark_matrices()