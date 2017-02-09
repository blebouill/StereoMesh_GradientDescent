#include <stdio.h>
#include <stdlib.h> 
#include <cmath>

#include "matrix_calculation_functions.h"
#include "other_functions_StereoMesh_GradientDescent.h"



// __ Inversion d'une matrice S (3x3)

int invert_3_3_matrix_S(double * S, double * inv_S)
{
	// Simplification de l'inversion car la matrice S contient des coordonnées homogènes 
	// donc la dernière ligne est composée de 1 (g = 1 ; h = 1; i = 1)

	double det_S;
	double a, b, c, d, e, f;
	double A, B, C, D, E, F, G, H, I;

	// Elements of S
	a = S[0];
	b = S[1];
	c = S[2];
	d = S[3];
	e = S[4];
	f = S[5];

	// Inverse elements
	A = e - f;
	B = f - d;
	C = d - e;
	D = c - b;
	E = a - c;
	F = b - a;
	G = b*f - c*e;
	H = c*d - a*f;
	I = a*e - b*d;

	// Determinant
	det_S = a*A + b*B + c*C;

	// Inversibility test
	if (det_S == 0.0)
	{
		printf("__ Matrice non inversible\n");
		return 1;
	}

	// Inverse matrix
	inv_S[0] = A / det_S;
	inv_S[1] = D / det_S;
	inv_S[2] = G / det_S;
	inv_S[3] = B / det_S;
	inv_S[4] = E / det_S;
	inv_S[5] = H / det_S;
	inv_S[6] = C / det_S;
	inv_S[7] = F / det_S;
	inv_S[8] = I / det_S;

	return 0;
}


// __ Inversion de toutes les matrices S

int invert_S_matrices(double * S, int N_T, double * inv_S)
{
	int i;
	double *ptr_S = NULL, *ptr_inv_S = NULL;

	// Récupération des pointeurs
	ptr_S = S;
	ptr_inv_S = inv_S;

	// Inversions
	for (i = 0; i < N_T; i++)
	{
		invert_3_3_matrix_S(ptr_S, ptr_inv_S);
		ptr_S += 9;
		ptr_inv_S += 9;
	}

	return 0;
}





// __ Calculs des S^(-1) x pH

int compute_inv_S_by_pH(int * img_label, int width_img, int height_img, double * inv_S, int N_T, double * inv_S_by_pH)
{
	int x, y, m, ind_pix = 0;
	double pH[3] = { 0.0, 0.0, 1.0 };

	double * ptr_out = inv_S_by_pH;


	for (y = 0; y < height_img; y++)
		for (x = 0; x < width_img; x++)
		{
			m = img_label[ind_pix];

			pH[0] = (double)x;
			pH[1] = (double)y;

			mult_3_3_matrix_by_3_1_vector(inv_S + 9*m, pH, ptr_out);

			ptr_out += 3;
			ind_pix++;
		}

	return 0;
}





// __ Calculs des matrices B = S^(-1) x K
int compute_B_matrices(double * inv_S, int N_T, double * K, double * B)
{
	int i;
	double * ptr_inv_S = inv_S, * ptr_B = B;

	for (i = 0; i < N_T; i++)
	{
		mult_two_3_3_matrices(ptr_inv_S, K, ptr_B);

		ptr_inv_S += 9;
		ptr_B += 9;
	}

	return 0;
}


// __ Calcul des matrices A = B x B'

int compute_A_matrices(double * B, int N_T, double * A)
{
	int i;
	double * ptr_A = A, * ptr_B = B;

	for (i = 0; i < N_T; i++)
	{
		mult_3_3_matrix_by_its_transpose(ptr_B, ptr_A);

		ptr_A += 9;
		ptr_B += 9;
	}

	return 0;
}



// __ Calcul de compute_nb_triangles_using_vertex
int compute_nb_triangles_using_vertex(int * ind_triangles_using_vertex, int N_V, int * nb_triangles_using_vertex)
{
	int i, j;

	int * ptr_nb_triangles_using_vertex = nb_triangles_using_vertex;
	int * ptr_ind_triangles_using_vertex = ind_triangles_using_vertex;

	// Initialization
	for (i = 0; i < N_V; i++)
	{
		*ptr_nb_triangles_using_vertex = 0;
		ptr_nb_triangles_using_vertex++;
	}

	// Compute nb_triangles_using_vertex
	for (i = 0; i < N_V; i++)
	{
		for (j = 0; j < 6; j++)
		{
			if (*ptr_ind_triangles_using_vertex >= 0)
				(*ptr_nb_triangles_using_vertex)++;

			ptr_ind_triangles_using_vertex++;
		}
		ptr_nb_triangles_using_vertex++;
	}

	return 0;
}












// __ Calculs préliminaires
int run_StereoMesh_precomputations(double * S, int N_T, double * inv_S)
{

	// __ Inversion des matrices S
	invert_S_matrices(S, N_T, inv_S);
}