#ifndef SYMNMF_HEADER
#define SYMNMF_HEADER
double** symnmf(double** H_mat, double** W_mat, int N, int k);
double** norm(double** A_mat, double** D_mat, int N);
double** ddg(double** A_mat, int N);
double** sym(double** X_mat, int N, int cols);
void freeMat(double** mat, int N);
#endif