#pragma once


// Multiplication de vecteurs
double mult_1_3_vector_by_3_1_vector(double * a, double * b);
int mult_3_1_vector_by_1_3_vector(double * a, double * b, double * C);


// Multiplication de vecteurs/matrices
int mult_3_3_matrix_by_3_1_vector(double * A, double * b, double * c);
int mult_1_3_vector_by_3_3_matrix(double * a, double * B, double * c);
int mult_1_3_vector_by_3_3_transpose_matrix(double * a, double * B, double * c);


// Multiplication de matrices
int mult_two_3_3_matrices(double * A, double * B, double * C);
int mult_3_3_matrix_by_its_transpose(double * A, double * B);


// Inversion de matrices
int invert_3_3_matrix(double * M, double * inv_M);

