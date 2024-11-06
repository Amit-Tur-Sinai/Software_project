# include <stdio.h>
# include <stdlib.h>
# include <math.h>


double find_dist(double* x1, double* x2, int cols) {
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
                dist = find_dist(X_mat[i], X_mat[j], cols);
                A_mat[i][j] = exp(-0.5 * dist);
            }
        }
    }

    return A_mat;
}

double** ddg(double** X_mat, int N) {
    double **D_mat;
    int i;

    D_mat = (double **)malloc(N*sizeof(double *));
    if (D_mat == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for (i=0;i<N;i++) {

    }

}