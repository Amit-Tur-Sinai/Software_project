# include <stdio.h>
# include <stdlib.h>
# include <math.h>

const double EPS = 0.0001;
const int max_iter = 300;
const double beta = 0.5;

double findDist(double* x1, double* x2, int cols) {
    int i;
    double dist = 0.0;
    
    for (i=0; i<cols;i++) {
        dist += pow(x1[i]-x2[i], 2);
    }
    
    return dist;
}


void freeMat(double** mat, int N) {
    int i;

    for (i=0;i<N;i++) {
        free(mat[i]);
    }
    free(mat);
}

double** MatrixMultiply(double** mat1, double** mat2, int rows1, int cols1, int rows2, int cols2) {
    double** result;
    int i, j, k;
    double sum;

    result = (double **)malloc(rows1*sizeof(double *));
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    // init result matrix in size rows1 * cols2
    for (i=0;i<rows1;i++) {
        result[i] = (double *)malloc(cols2*sizeof(double));
    }

    for (i=0;i<rows1;i++) {
        for (j=0;j<cols2;j++) {
            sum = 0.0; // init sum
            for (k=0;k<rows2;k++) {
                sum += mat1[i][k] * mat2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

double squaredFrobNorm(double** mat, int N, int cols) {
    int i, j;
    double sum = 0.0;
    for (i=0;i<N;i++) {
        for (j=0;j<cols;j++) {
            sum += pow(mat[i][j], 2);
        }
    }
    return sum;
}

int covergence(double** H_old, double** H_new, int N, int cols) {
    int i, j, output;
    double** result;
    double frobNorm;
    result = (double **)malloc(N*sizeof(double *));
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<N;i++) {
        result[i] = (double *)malloc(N*sizeof(double));
        if (result[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    for (i=0;i<N;i++) {
        for (j=0;j<cols;j++) {
            result[i][j] = H_new[i][j] - H_old[i][j];
        }
    }
    frobNorm = squaredFrobNorm(result, N, cols);
    if (frobNorm < EPS) {output = 1;}
    else {output = 0;}
    freeMat(result, N);
    return 0;
}

double** transpose(double** mat, int N, int cols) {
    double** result;
    int i, j;
    result = (double **)malloc(N*sizeof(double *));
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<N;i++) {
        result[i] = (double *)malloc(N*sizeof(double));
        if (result[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    for (i=0;i<N;i++) {
        for (j=0;j<cols;j++) {
            result[i][j] = mat[j][i];
        }
    }
    return result;
}

double** updateH(double** H_mat, double** W_mat, int N, int cols) {
    double **result, **numerator, **denominator, **HT;
    int i, j;
    result = (double **)malloc(N*sizeof(double *));
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<N;i++) {
        result[i] = (double *)malloc(N*sizeof(double));
        if (result[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    numerator = MatrixMultiply(W_mat, H_mat, N, N, N, cols);
    HT = transpose(H_mat, N, cols);
    denominator = MatrixMultiply(H_mat, HT, N, cols, cols, N);
    denominator = MatrixMultiply(denominator, H_mat, N, N, N, cols);
    for (i=0;i<N;i++) {
        for (j=0;j<cols;j++) {
            result[i][j] = 1-beta + beta*(numerator[i][j]/denominator[i][j]);
            result[i][j] = result[i][j] * H_mat[i][j];
        }
    }
    freeMat(HT, cols);
    freeMat(numerator, N);
    freeMat(denominator, N);
    return result;
}

double** sym(double** X_mat, int N, int cols) {
    double **A_mat;
    int i, j, dist;

    A_mat = (double **)malloc(N*sizeof(double *));
    if (A_mat == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for (i=0;i<N;i++) {
        A_mat[i] = (double *)malloc(N*sizeof(double));
        if (A_mat[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }

        for(j=0;j<N;j++) {
            if (i == j) {
                A_mat[i][j] = 0.0;
            }
            else {
                dist = findDist(X_mat[i], X_mat[j], cols);
                A_mat[i][j] = exp(-0.5 * dist);
            }
        }
    }
    return A_mat;
}

double** ddg(double** A_mat, int N) {
    double **D_mat;
    int i, j;
    double sum;
    
    D_mat = (double **)malloc(N*sizeof(double *));
    if (D_mat == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for (i=0;i<N;i++) {
        D_mat[i] = (double *)calloc(N,sizeof(double));
        if (D_mat[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
        sum = 0.0; // init sum variable
        for(j=0;j<N;j++) {
            sum += A_mat[i][j];
        }
        D_mat[i][i] = sum;
    }

    return D_mat;
}

double** norm(double** A_mat, double** D_mat, int N) {
    int i;
    double **tempMat, **result;

    // convert D_mat to be D_mat ** -0.5
    for(i=0;i<N;i++) {
        D_mat[i][i] = pow(D_mat[i][i], -0.5);
    }

    tempMat = MatrixMultiply(D_mat, A_mat, N, N, N, N);
    result = MatrixMultiply(tempMat, D_mat, N, N, N, N);

    free(tempMat);
    return result;
}

double** symnmf(double** H_mat, double** W_mat, int N, int k) {
    // iter 300 times, calling update and convergence
    updateH(H_mat, W_mat, N, k)
}


