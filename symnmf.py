import numpy as np
import pandas as pd
import sys
import symnmf 

np.random.seed(1234)

def init_H_mat(W_mat, N, k):
    m = np.mean(W_mat)
    high = 2* np.sqrt(m/k)
    H_mat = np.random.uniform(0, high, size=(N,k))
    return H_mat

if __name__ == '__main__':
    if len(sys.argv) == 4:
        k, goal, file_name = sys.argv[1:]

        X_mat = pd.read_csv(file_name, sep=",", header=None).values.tolist()
        N = len(X_mat)
        cols = len(X_mat[0])

        if goal == "symnmf":
            W_mat = symnmf.norm(X_mat)
            H_mat = init_H_mat(W_mat, N, k)
            result_mat = symnmf.symnmf(H_mat, W_mat)

        elif goal == "sym":
            A_mat = symnmf.sym(X_mat)

        elif goal == "ddg":
            D_mat = symnmf.ddg(X_mat)

        elif goal == "norm":
            W_mat = symnmf.norm(X_mat)

        else:
            print("An Error Has Occurred")
            exit()    
    else:
        print("An Error Has Occurred")
        exit()