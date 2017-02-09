#pragma once


// Inversion des matrices S
int invert_3_3_matrix_S(double * S, double * inv_S);
int invert_S_matrices(double * S, int N_T, double * inv_S);


// Calculs des S^(-1)*pH
int compute_inv_S_by_pH(int * img_label, int width_img, int height_img, double * inv_S, int N_T, double * inv_S_by_pH);


// Calcul des matrices B
int compute_B_matrices(double * inv_S, int N_T, double * K, double * B);


// Calcul des matrices A
int compute_A_matrices(double * B, int N_T, double * A);


// Calcul de compute_nb_triangles_using_vertex
int compute_nb_triangles_using_vertex(int * ind_triangles_using_vertex, int N_V, int * nb_triangles_using_vertex);
