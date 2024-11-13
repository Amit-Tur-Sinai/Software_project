# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

const double EPS = 0.0001;
const int MAX_ITER = 300;
const double BETA = 0.5;

double findDist(double* x1, double* x2, int cols) {
    int i;
    double dist = 0.0;
    
    for (i=0; i<cols;i++) {
        dist += pow(x1[i]-x2[i], 2);
    }
    
    return dist;
}

/*Free matrix memory*/
void freeMat(double** mat, int N) {
    int i;

    for (i=0;i<N;i++) {
        free(mat[i]);
    }
    free(mat);
}

void printMat(double** mat, int N, int cols) {
    int i, j;
    for (i=0;i<N;i++) {
        for (j=0;j<cols;j++) {
            printf("%.4f", mat[i][j]);
            if (j != cols -1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

/* Multiplying 2 matrices */
double** MatrixMultiply(double** mat1, double** mat2, int rows1, int cols1, int rows2, int cols2) {
    double** result;
    int i, j, k;
    double sum;

    if (cols1 != rows2) { /*Matrices dimensions do not align */
        printf("An Error Has Occurred\n");
        exit(1);
    }
    /*Initialize result matrix*/
    result = (double **)malloc(rows1*sizeof(double *));
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<rows1;i++) {
        result[i] = (double *)malloc(cols2*sizeof(double));
        if (result[i] == NULL) {
            printf("An Error Has Occurred\n");
            for (j=0;j<i;j++) {
                free(result[j]);
            }
            free(result);
            exit(1);
        }
    }
    for (i=0;i<rows1;i++) {
        for (j=0;j<cols2;j++) {
            sum = 0.0; /* init sum */
            for (k=0;k<rows2;k++) {
                sum += mat1[i][k] * mat2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

/* Implementation of the Frobenius norm */
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

/* Checking convergence of two matrices by using the Frobenius norm*/
int convergence(double** H_old, double** H_new, int N, int cols) {
    int i, j, output;
    double** result;
    double frobNorm;
    /*Initialize result matrix*/
    result = (double **)malloc(N*sizeof(double *));
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<N;i++) {
        result[i] = (double *)malloc(cols*sizeof(double));
        if (result[i] == NULL) {
            printf("An Error Has Occurred\n");
            for (j=0;j<i;j++) {
                free(result[j]);
            }
            free(result);
            exit(1);
        }
    }
    for (i=0;i<N;i++) { /*Calculating the difference matrix*/
        for (j=0;j<cols;j++) {
            result[i][j] = H_new[i][j] - H_old[i][j];
        }
    }
    frobNorm = squaredFrobNorm(result, N, cols); /*Convergence check of the difference matrix*/
    if (frobNorm < EPS) {
        output = 1;
    }
    else {
        output = 0;
    }
    freeMat(result, N);
    return output;
}

/* Transposing a matrix by allocating space for a new matrix, iterating over the original and changing the indexes */
double** transpose(double** mat, int N, int cols) {
    int i, j;
    
    /*Initialize result matrix*/
    double** result = (double **)malloc(cols * sizeof(double *));
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for (i = 0; i < cols; i++) {
        result[i] = (double *)malloc(N * sizeof(double));
        if (result[i] == NULL) {
            printf("An Error Has Occurred\n");
            /* Free previously allocated rows */
            for (j=0;j<i;j++) {
                free(result[j]);
            }
            free(result);
            exit(1);
        }
    }
    
    for (i=0;i<N;i++) {
        for (j=0;j<cols;j++) {
            result[j][i] = mat[i][j];
        }
    }

    return result;
}

/* Updating H matrix according to the recursive formula */
double** updateH(double** H_mat, double** W_mat, int N, int cols) {
    double **result, **numerator, **denominator, **HT;
    int i, j;
    /*Initialize result matrix*/
    result = (double **)malloc(N*sizeof(double *));
    if (result == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<N;i++) {
        result[i] = (double *)malloc(cols*sizeof(double));
        if (result[i] == NULL) {
            printf("An Error Has Occurred\n");
            for (j=0;j<i;j++) {
                free(result[j]);
            }
            free(result);
            exit(1);
        }
    }
    /* Get the numerator and denominator of the expression by multiplying matrices, and transposing */
    numerator = MatrixMultiply(W_mat, H_mat, N, N, N, cols);
    HT = transpose(H_mat, N, cols);
    denominator = MatrixMultiply(H_mat, HT, N, cols, cols, N);
    denominator = MatrixMultiply(denominator, H_mat, N, N, N, cols);
    /* Using the numerator and denominator matrices, assign the values according to the formula */
    for (i=0;i<N;i++) {
        for (j=0;j<cols;j++) {
            result[i][j] = 1-BETA + BETA*(numerator[i][j]/denominator[i][j]);
            result[i][j] = result[i][j] * H_mat[i][j];
        }
    }
    freeMat(HT, cols);
    freeMat(numerator, N);
    freeMat(denominator, N);
    return result;
}

/* SYM Implementation: Using the findDist helper function, we calculate the distance of every two points, and assign the value to the newly allocated matrix */
double** sym(double** X_mat, int N, int cols) {
    double **A_mat;
    int i, j;
    double dist;

    A_mat = (double **)malloc(N*sizeof(double *));
    if (A_mat == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<N;i++) {
        A_mat[i] = (double *)malloc(N*sizeof(double));
        if (A_mat[i] == NULL) {
            printf("An Error Has Occurred\n");
            for (j=0;j<i;j++) {
                free(A_mat[j]);
            }
            free(A_mat);
            exit(1);
        }
    }
    /* Calculating the A matrix*/
    for (i=0;i<N;i++) {
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

/* DDG implementation: Using the similarity matrix, we allocate space for a new matrix, calculate the values of the lines of A, and assign the values to the new matrix */
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
            for (j=0;j<i;j++) {
                free(D_mat[j]);
            }
            free(D_mat);
            exit(1);
        }
    }
    /* Calculating the D matrix*/
    for (i=0;i<N;i++) {        
        sum = 0.0; /* init sum variable */
        for(j=0;j<N;j++) {
            sum += A_mat[i][j];
        }
        D_mat[i][i] = sum;
    }

    return D_mat;
}

/* NORM implementation: Using the multiply helper function, we calculate the norm matrix */
double** norm(double** A_mat, double** D_mat, int N) {
    int i;
    double **tempMat, **result;

    /* convert D_mat to be D_mat ** -0.5 */
    for(i=0;i<N;i++) {
        D_mat[i][i] = pow(D_mat[i][i], -0.5);
    }

    tempMat = MatrixMultiply(D_mat, A_mat, N, N, N, N);
    result = MatrixMultiply(tempMat, D_mat, N, N, N, N);

    freeMat(tempMat, N);
    return result;
}

/* SYMNMF implementation: Iterate up to 300 times, and check convergence, using the two helper functions updateH() and convergence() */
double** symnmf(double** H_mat, double** W_mat, int N, int k) {
    int iter, end;
    double** H_temp;
    double** H_new = H_mat;

    for (iter=0;iter<MAX_ITER;iter++) {
        H_temp = H_new;  
        H_new = updateH(H_temp, W_mat, N, k);  
        end = convergence(H_temp, H_new, N, k); 

        if (H_temp != H_mat) { /* Free H_temp if it is not the original input matrix */       
            freeMat(H_temp, N);
        }

        H_mat = H_new; /* Updates H_mat pointer*/
        if (end) { /* Break the loop if convergence is achieved */
            break;
        }
    }
    return H_new;
}

int main(int argc, char *argv[]) {
    FILE *file;
    int i, j, N, cols;
    char *line, *val, *goal, *fileName; 
    double **X_mat, **A_mat, **W_mat, **D_mat;
    size_t len;

    if (argc != 3) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    goal = argv[1];
    fileName = argv[2];
    file = fopen(fileName, "r");
    if (file == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    N = 0;
    cols = 0;
    
    /* Calculating the dimensions of the data */
    line = NULL;
    len = 0;
    while (getline(&line, &len, file) != -1) {
        if (N == 0) {
            val = strtok(line, ",");
            while (val != NULL) {
                cols++;
                val = strtok(NULL, ",");
            }
        }
        N++;
    }
    rewind(file); /* Going back to the beginning of the file in order to read the points */

    /*Initialize X matrix*/
    X_mat = (double **)malloc(N*sizeof(double *));
    if (X_mat == NULL) {
        printf("An Error Has Occurred\n");
        free(line);
        fclose(file);
        exit(1);
    }
    for (i=0;i<N;i++) {
        X_mat[i] = (double *)malloc(N*sizeof(double));
        if (X_mat[i] == NULL) {
            printf("An Error Has Occurred\n");
            for (j=0;j<i;j++) {
                free(X_mat[j]);
            }
            free(X_mat);
            free(line);
            fclose(file);
            exit(1);
        }
    }
    /*Reading the data from the file to the X matrix*/
    for (i=0;i<N;i++) {
        if (getline(&line, &len, file) != -1) {
            val = strtok(line, ",");
            for (j=0;j<cols;j++) {
                if (val != NULL) {
                    X_mat[i][j] = atof(val);
                    val = strtok(NULL, ",");
                }
            }
        }
    }
    free(line);
    fclose(file);

    /* Start of logic */
    if (strcmp(goal, "sym") == 0) {
        A_mat = sym(X_mat, N, cols);
        printMat(A_mat, N, N);
        freeMat(X_mat, N);
        freeMat(A_mat, N);
    }
    else if (strcmp(goal, "ddg") == 0) {
        A_mat = sym(X_mat, N, cols);
        D_mat = ddg(A_mat, N);
        printMat(D_mat, N, N);
        freeMat(X_mat, N);
        freeMat(A_mat, N);
        freeMat(D_mat, N);
    }

    else if (strcmp(goal, "norm") == 0) {
        A_mat = sym(X_mat, N, cols);
        D_mat = ddg(A_mat, N);
        W_mat = norm(A_mat, D_mat, N);
        printMat(W_mat, N, N);
        freeMat(X_mat, N);
        freeMat(A_mat, N);
        freeMat(D_mat, N);
        freeMat(W_mat, N);
    }
    else {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    return 0;
}