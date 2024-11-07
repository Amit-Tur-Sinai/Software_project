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

double** MatrixMultiply(double** mat1, double** mat2, int rows1, int cols1, int rows2, int cols2) {
    // Multiplying matrices
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
    // Implementation of the Frobenius norm
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
    // Checking convergence of two matrices by using the Frobenius norm
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
    // Transposing a matrix by allocating space for a new matrix, iterating over the original and changing the indexes
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
    // Helper function, we use it in symnmf
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
    // Get the numerator and denominator of the expression by multiplying matrices, and transposing
    numerator = MatrixMultiply(W_mat, H_mat, N, N, N, cols);
    HT = transpose(H_mat, N, cols);
    denominator = MatrixMultiply(H_mat, HT, N, cols, cols, N);
    denominator = MatrixMultiply(denominator, H_mat, N, N, N, cols);
    // Using the numerator and denominator matrices, assign the values according to the formula
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

double** sym(double** X_mat, int N, int cols) {
    // Using the findDist helper function, we calculate the distance of every two points, and assign the value to the newly allocated matrix
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
    // Using the similarity matrix, we allocate space for a new matrix, calculate the values of the lines of A, and assign the values to the new matrix
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
    // Using the transpose and multiply helper functions, we calculate the norm matrix
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
    // Iterate up to 300 times, and check convergence, using the two helper functions
    int iter, i, end;
    double** H_new;
    H_new = (double **)malloc(N*sizeof(double *));
    if (H_new == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<N;i++) {
        H_new[i] = (double *)malloc(N*sizeof(double));
        if (H_new[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    for (iter=0;iter<MAX_ITER;iter++) {
        H_new = updateH(H_mat, W_mat, N, k);
        end = covergence(H_mat, H_new, N, k);
        H_mat = H_new;
        if (end) {
            break;
        }
    }
    freeMat(H_new, N);
    return H_mat;
}

int main(int argc, char *argv[]) {
    FILE *file;
    int i, j, N, cols;
    char *line, *val; 
    double **X_mat, **A_mat, **W_mat, **D_mat;
    if (argc != 3) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    char* goal = argv[1];
    char* fileName = argv[2];
    file = fopen(fileName, "r");
    if (file == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    N = 0;
    cols = 0;
    
    // Calculating the dimensions of the data
    size_t len = 0;
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
    rewind(file); // Going back to the beginning of the file in order to read the points

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

    // Start of logic
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


