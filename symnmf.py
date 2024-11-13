import numpy as np
import pandas as pd
import sys
import symnmf 

# Sets numpy's random seed
np.random.seed(1234)

# Initialize H matrix using the W matrix
def init_H_mat(W_mat, N, k):
    m = np.mean(W_mat)
    high = 2* np.sqrt(m/k)
    H_mat = np.random.uniform(0, high, size=(N,k))
    return H_mat.tolist()

def printResult(resultMat):
    for row in resultMat:
        print(','.join([f'{num:.4f}' for num in row]))  

if __name__ == '__main__':
    try:
        if len(sys.argv) == 4:
            # Getting the arguments from the CMD
            k, goal, file_name = sys.argv[1:]

            # Converting the text file to list using the pandas package
            X_mat = pd.read_csv(file_name, sep=",", header=None).values.tolist()
            N = len(X_mat)
            cols = len(X_mat[0])
            k = int(k)

            if goal == "symnmf":
                W_mat = symnmf.norm(X_mat)
                H_mat = init_H_mat(W_mat, N, k)
                result_mat = symnmf.symnmf(H_mat, W_mat, N, k)
                printResult(result_mat)    

            elif goal == "sym":
                A_mat = symnmf.sym(X_mat)
                printResult(A_mat)                

            elif goal == "ddg":
                D_mat = symnmf.ddg(X_mat)
                printResult(D_mat)                

            elif goal == "norm":
                W_mat = symnmf.norm(X_mat)
                printResult(W_mat)    

            else:
                print("An Error Has Occurred")
                exit() 
        # Else the number of arguments is not 4 as expected    
        else:
            print("An Error Has Occurred")
            exit()
    except:
        print("An Error Has Occurred")
        exit()         