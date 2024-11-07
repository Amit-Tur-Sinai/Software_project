# include <stdio.h>
# include <stdlib.h>
# include <math.h>



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
    
    // init reault matrix in size rows1 * cols2
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


