import qldpc
import galois
import numpy as np
import json

def load_boundaries(filename):
    print("checking {:} for boundary operators.".format(filename))
    with open(filename, 'r') as f:
        gf = galois.GF2
        raw_data = json.load(f)
        print("n_rows x n_cols = {:} x {:} ".format(raw_data["boundary_up_n_rows"], raw_data["boundary_up_n_cols"]))
        boundary_up_raw = raw_data["boundary_up"]
        boundary_up = []
        ix = 0
        for row in range(raw_data["boundary_up_n_rows"]):
            new_row = []
            for col in range(raw_data["boundary_up_n_cols"]):
                new_row.append(boundary_up_raw[ix])
                ix += 1
            boundary_up.append(new_row)
        boundary_down_raw = raw_data["boundary_down"]
        boundary_down = []
        ix = 0
        for row in range(raw_data["boundary_down_n_rows"]):
            new_row = []
            for col in range(raw_data["boundary_down_n_cols"]):
                new_row.append(boundary_down_raw[ix])
                ix += 1
            boundary_down.append(new_row)
    return (gf(boundary_down), gf(boundary_up))

        
if __name__ == "__main__":
    boundary_down, boundary_up = load_boundaries('/Users/matt/repos/qec/tmp/boundary_operators.json')
    
    boundary_down = boundary_down.T
    boundary_up = boundary_up.T
    print('boundary down shape: ', boundary_down.shape)
    print('boundary up shape: ', boundary_up.shape)
    
    good_boundary_down_rows = []
    for row in range(boundary_down.shape[0]):
        row_product = boundary_down[row,:] @ boundary_up.T
        if not np.any(row_product):
            good_boundary_down_rows.append(row)
        # print(row_product)
    print("good rows: ", good_boundary_down_rows)
    boundary_down = boundary_down[good_boundary_down_rows,:]
    quantum_code = qldpc.codes.CSSCode(boundary_down, boundary_up)
    print("Code constructed!")
    print(quantum_code.get_code_params())