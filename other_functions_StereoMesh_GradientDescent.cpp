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



// __ Calcul de la carte des disparités
int compute_disparity_map(int * img_label, int N_pixels, double * D, double * inv_S_by_pH, double * disparity_map)
{
	int i, m;

	int * ptr_img_label = img_label;
	double * ptr_inv_S_by_pH = inv_S_by_pH;
	double * ptr_disparity_map = disparity_map;

	for (i = 0; i < N_pixels; i++)
	{
		m = *ptr_img_label;
		*ptr_disparity_map = mult_1_3_vector_by_3_1_vector(D + 3 * m, ptr_inv_S_by_pH);		// disp_pix = d_m x S_m^(-1) x pH

		ptr_img_label++;
		ptr_inv_S_by_pH++;
		ptr_disparity_map++;
	}

	return 0;
}



// __ Calcul de l'image interpolée (RGB)
int compute_img2_interp(double * img1, int width_img, int height_img, double * disparity_map, double * img2_interp)
{
	// Ne prend pas en compte si des pixels sont interpolés plusieurs fois

	int i, x, y, x1, N_pixels = height_img * width_img;
	double x_img2;

	double * ptr_disparity_map = disparity_map;

	double * ptr_img1_R = img1;
	double * ptr_img1_G = img1 + N_pixels;
	double * ptr_img1_B = img1 + 2 * N_pixels;

	double * ptr_img2_interp = img2_interp;
	double * ptr_img2_interp_R = img2_interp;
	double * ptr_img2_interp_G = img2_interp + N_pixels;
	double * ptr_img2_interp_B = img2_interp + 2 * N_pixels;

	// Initialization
	for (i = 0; i < 3 * N_pixels; i++)
	{
		*ptr_img2_interp = 0.0;
		ptr_img2_interp++;
	}
		
	// Loop
	for (y = 0; y < height_img; y++)
	{
		for (x = 0; x < width_img; x++)
		{
			// Pixel disparity
			x_img2 = (double)x + *ptr_disparity_map;

			if ((x_img2 >= 0.0) && (x_img2 <= (double)(width_img - 1)))
			{
				// Interpolations (neirest neighbor)
				x1 = (int)round(x_img2);

				ptr_img2_interp_R[x1] = *ptr_img1_R;
				ptr_img2_interp_G[x1] = *ptr_img1_G;
				ptr_img2_interp_B[x1] = *ptr_img1_B;
			}

			// Update pointers (next pixel)
			ptr_disparity_map++;
			ptr_img1_R++;
			ptr_img1_G++;
			ptr_img1_B++;
		}

		// Update pointers (next line)
		ptr_img2_interp_R += width_img;
		ptr_img2_interp_G += width_img;
		ptr_img2_interp_B += width_img;


	}

	return 0;
}



// __ Calcul du nouveau gradient
int compute_next_grad_D(double * grad_D_DATA, double * grad_D_BREACH, double * grad_D_NORMAL, int N_T, double lambda_BREACH, double lambda_NORMAL, double * grad_D_next)
{
	int i;

	double * ptr_grad_D_DATA = grad_D_DATA;
	double * ptr_grad_D_BREACH = grad_D_BREACH;
	double * ptr_grad_D_NORMAL = grad_D_NORMAL;
	double * ptr_grad_D_next = grad_D_next;

	for (i = 0; i < 3 * N_T; i++)
	{
		*ptr_grad_D_next = (*ptr_grad_D_DATA) + lambda_BREACH * (*ptr_grad_D_BREACH) + lambda_NORMAL * (*ptr_grad_D_NORMAL);

		ptr_grad_D_DATA++;
		ptr_grad_D_BREACH++;
		ptr_grad_D_NORMAL++;
		ptr_grad_D_next++;

	}

	return 0;
}



// __ Mise à jour de D
int update_D(double * D, double * grad_D, int N_T, double delta_GradientDescent, double * D_new)
{
	int i;

	double * ptr_D = D;
	double * ptr_grad_D = grad_D;
	double * ptr_D_new = D_new;

	for (i = 0; i < 3 * N_T; i++)
	{
		*ptr_D_new = (*ptr_D) - delta_GradientDescent * (*ptr_grad_D);

		ptr_D_new++;
		ptr_D++;
		ptr_grad_D++;
	}

	return 0;
}