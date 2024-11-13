import math
import symnmf 
import sys
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score

EPS = 0.0001
MAX_ITER = 300

np.random.seed(1234)

# ######################### HW1 KMeans Part #########################

def findDist(t1, t2):
    s = 0
    for i in range(len(t1)):
        s += (t1[i] - t2[i]) ** 2
    return math.sqrt(s)

def check_convergence(lst):
    for centroid in lst:
        if not centroid[2]: # If there is at least one centroid that changed in the last iteration then we want to go again
            return False
    return True

def update_centroid(centroid):
    prev_centroid = centroid[0].copy()
    for i in range(centroid[3]):
        s = 0
        for member in centroid[1]:
            s += member[i]
        centroid[0][i] = s / len(centroid[1]) # The new centroid value

    if len(prev_centroid) != 0:
        centroid[2] = (findDist(centroid[0], prev_centroid) < EPS)

    centroid[1] = [] # reset the members list

def kmeans(k, iter, dataArray):   
    centroids = [] 
    for i in range(k): # Setting the first k data points as the initial cenotroids
        # centroid value, list of the members of its class, boolean indicator for convergence, point dimension
        centroid = [dataArray[i].copy(), [], False, len(dataArray[i])] 
        centroids.append(centroid)

    for i in range(iter):
        if check_convergence(centroids):
            break

        for point in dataArray:
            min_dist = float('inf')
            candidate_centroid = None
            for centroid in centroids:
                dist = findDist(centroid[0], point)
                if dist < min_dist:
                    candidate_centroid = centroid
                    min_dist = dist
            candidate_centroid[1].append(point)

        for centroid in centroids:
            update_centroid(centroid)
    
    return centroids

# ######################### Silhouette Analysis Part #########################

# Initialize H matrix using the W matrix
def init_H_mat(W_mat, N, k):
    m = np.mean(W_mat)
    high = 2* np.sqrt(m/k)
    H_mat = np.random.uniform(0, high, size=(N,k))
    return H_mat.tolist()

# Generates symnmf clusters using the symnmf module
def generate_symnmf_clusters(X_mat, k):
    W_mat = symnmf.norm(X_mat)
    H_mat = init_H_mat(W_mat, N, k)
    result_mat = symnmf.symnmf(H_mat, W_mat, N, k)
    
    clusters = np.argmax(result_mat, axis=1) # choose for each element the cluster with the highest association score.
    return clusters

# Generates kmeans clusters using the kmeans implemantation from HW1
def generate_kmeans_clusters(X_mat, k):
    centroids = kmeans(k, MAX_ITER, X_mat)
    clusters = []

    # For each data point, choose the nearest centroid to create clusters
    for point in X_mat:
        min_dist = float('inf')
        candidate_idx = None
        for idx, centroid in enumerate(centroids): # Iterating over all centroids
            dist = findDist(centroid[0], point)
            if dist < min_dist:
                candidate_idx = idx
                min_dist = dist
        clusters.append(candidate_idx)
    
    return clusters

if __name__ == '__main__':
    try:
        if len(sys.argv) == 3:
            # Getting the arguments from the CMD
            k, file_name = sys.argv[1:]

            # Converting the text file to list using the pandas package
            X_mat = pd.read_csv(file_name, sep=",", header=None).values.tolist()
            N = len(X_mat)
            cols = len(X_mat[0])
            k = int(k)

            symnmf_clusters = generate_symnmf_clusters(X_mat, k)
            kmeans_clusters = generate_kmeans_clusters(X_mat, k)
            
            symnmf_silhouette_score = silhouette_score(X_mat, symnmf_clusters)
            kmeans_silhouette_score = silhouette_score(X_mat, kmeans_clusters)

            print(f"nmf: {symnmf_silhouette_score:.4f}")
            print(f"kmeans: {kmeans_silhouette_score:.4f}")
        
        # Else the number of arguments is not 3 as expected  
        else:
            print("An Error Has Occurred")
            exit()             
    except:
        print("An Error Has Occurred")
        exit()   
