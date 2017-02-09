#include <stdio.h>
#include <stdlib.h> 
#include <cmath>

#include "matrix_calculation_functions.h"
#include "energy_StereoMesh_GradientDescent.h"


//===========================================================================
///	compute_grad_E_DATA
///
/// Compute gradient and current energy of the E_DATA term
//===========================================================================
int compute_grad_E_DATA(
			double *	img1,					// Image 1	(RGB)
			double *	img2,					// Image 2	(RGB)
			double *	grad_img2_x,			// Image 2 gradient along x (RGB)
			int			width_img,				// Image width
			int			height_img,				// Image height
			int			N_pixels,				// N_pixels = height_img x width_img;
			int *		img_label,				// Triangulation labels (segmentation result) [0 to N_T-1]
			int			N_T,					// Number of triangles
			double *	inv_S_by_pH,			// S^(-1) x pH
			double *	D,						// Disparity vector
			double *	grad_D_DATA,			// Gradient of the E_DATA term
			double *	current_energy_E_DATA	// Current energy of the E_DATA term
)
{
	// =========================================
	//	Variables

	int i, x, y, m, x1, x2, ind_pix_x1;
	double x_img2, alpha;
	double delta21_grad_R, delta21_grad_G, delta21_grad_B, sum_delta;
	double delta_img21_R, delta_img21_G, delta_img21_B;
	double interp_img2_R, interp_img2_G, interp_img2_B;
	double interp_grad_img2_x_R, interp_grad_img2_x_G, interp_grad_img2_x_B;

	int * norm_grad_nb_pix = NULL;
	int * ptr_norm_grad_nb_pix = NULL;

	int * ptr_img_label = img_label;

	double * ptr_inv_S_by_pH = inv_S_by_pH;
	double * ptr_grad_D_DATA = grad_D_DATA;

	double * ptr_img1_R = img1;
	double * ptr_img1_G = img1 + N_pixels;
	double * ptr_img1_B = img1 + 2 * N_pixels;

	double * ptr_img2_R = img2;
	double * ptr_img2_G = img2 + N_pixels;
	double * ptr_img2_B = img2 + 2 * N_pixels;

	double * ptr_grad_img2_x_R = grad_img2_x;
	double * ptr_grad_img2_x_G = grad_img2_x + N_pixels;
	double * ptr_grad_img2_x_B = grad_img2_x + 2 * N_pixels;



	// =========================================
	//	Memory Allocation

	norm_grad_nb_pix = (int *)calloc(N_T, sizeof(int));



	// =========================================
	//	Main Loop

	// Initialization
	for (i = 0; i < 3 * N_T; i++)
	{
		*ptr_grad_D_DATA = 0.0;
		ptr_grad_D_DATA++;
	}

	*current_energy_E_DATA = 0.0;

	// Loop
	for (y = 0; y < height_img; y++)
	{
		for (x = 0; x < width_img; x++)
		{
			// Pixel label
			m = *ptr_img_label;

			// Pixel disparity
			x_img2 = (double)x + mult_1_3_vector_by_3_1_vector(D + 3*m, ptr_inv_S_by_pH);

			if ((x_img2 >= 0.0) && (x_img2 <= (double)(width_img - 1)))
			{
				// Interpolations (linear)
				alpha = x_img2 - floor(x_img2);
				x1 = (int)floor(x_img2);
				x2 = (int)floor(x_img2 + 1);

				if (x2 < width_img)
				{
					interp_img2_R = alpha * ptr_img2_R[x1] + (1 - alpha) * ptr_img2_R[x2];
					interp_img2_G = alpha * ptr_img2_G[x1] + (1 - alpha) * ptr_img2_G[x2];
					interp_img2_B = alpha * ptr_img2_B[x1] + (1 - alpha) * ptr_img2_B[x2];

					interp_grad_img2_x_R = alpha * ptr_grad_img2_x_R[x1] + (1 - alpha) * ptr_grad_img2_x_R[x2];
					interp_grad_img2_x_G = alpha * ptr_grad_img2_x_G[x1] + (1 - alpha) * ptr_grad_img2_x_G[x2];
					interp_grad_img2_x_B = alpha * ptr_grad_img2_x_B[x1] + (1 - alpha) * ptr_grad_img2_x_B[x2];
				}
				else	// x2 > (width_img - 1) , i.e. x1 = (width_img - 1)
				{
					interp_img2_R = ptr_img2_R[x1];
					interp_img2_G = ptr_img2_G[x1];
					interp_img2_B = ptr_img2_B[x1];

					interp_grad_img2_x_R = ptr_grad_img2_x_R[x1];
					interp_grad_img2_x_G = ptr_grad_img2_x_G[x1];
					interp_grad_img2_x_B = ptr_grad_img2_x_B[x1];
				}

				// Gradient
				delta_img21_R = (interp_img2_R - *ptr_img1_R);
				delta_img21_G = (interp_img2_G - *ptr_img1_G);
				delta_img21_B = (interp_img2_B - *ptr_img1_B);

				delta21_grad_R = 2.0 * delta_img21_R * interp_grad_img2_x_R;
				delta21_grad_G = 2.0 * delta_img21_G * interp_grad_img2_x_G;
				delta21_grad_B = 2.0 * delta_img21_B * interp_grad_img2_x_B;

				sum_delta = delta21_grad_R + delta21_grad_G + delta21_grad_B;

				grad_D_DATA[3 * m]     += sum_delta * ptr_inv_S_by_pH[0];
				grad_D_DATA[3 * m + 1] += sum_delta * ptr_inv_S_by_pH[1];
				grad_D_DATA[3 * m + 2] += sum_delta * ptr_inv_S_by_pH[2];

				// Update energy
				*current_energy_E_DATA += delta_img21_R*delta_img21_R + delta_img21_G*delta_img21_G + delta_img21_B*delta_img21_B;

				// Normalization vector
				norm_grad_nb_pix[m]++;
			}

			// Update pointers (next pixel)
			ptr_img_label++;

			ptr_inv_S_by_pH += 3;

			ptr_img1_R++;
			ptr_img1_G++;
			ptr_img1_B++;
		}

		// Update pointers (next line)
		ptr_img2_R += width_img;
		ptr_img2_G += width_img;
		ptr_img2_B += width_img;

		ptr_grad_img2_x_R += width_img;
		ptr_grad_img2_x_G += width_img;
		ptr_grad_img2_x_B += width_img;

	}


	
	// =========================================
	//	Normalization

	ptr_grad_D_DATA = grad_D_DATA;
	ptr_norm_grad_nb_pix = norm_grad_nb_pix;

	for (i = 0; i < N_T; i++)
	{
		if (*ptr_norm_grad_nb_pix != 0)
		{
			ptr_grad_D_DATA[0] /= 3.0 * (double)(*ptr_norm_grad_nb_pix);	// 3.0 for the number of color channels
			ptr_grad_D_DATA[1] /= 3.0 * (double)(*ptr_norm_grad_nb_pix);
			ptr_grad_D_DATA[2] /= 3.0 * (double)(*ptr_norm_grad_nb_pix);
		}
		ptr_norm_grad_nb_pix++;
		ptr_grad_D_DATA += 3;
	}



	// =========================================
	//	Memory Deallocation

	free(norm_grad_nb_pix);

	return 0;
}










//===========================================================================
///	compute_grad_E_BREACH
///
/// Compute gradient and current energy of the E_BREACH term
//===========================================================================
int compute_grad_E_BREACH(
			int *		ind_vertex_in_triangle,		// Indices of the 3 vertices in each triangle [size 3*N_T]
			int *		ind_triangles_using_vertex,	// Indices of the triangles using each vertex [size 6*N_V] (max number of triangles using a vertex is 6, when less the last indices are -1)
			int *		nb_triangles_using_vertex,	// Number of triangles using each vertex [size N_V]
			int			N_V,						// Number of vertices
			int			N_T,						// Number of triangles
			double *	D,							// Disparity vector
			double *	grad_D_BREACH,				// Gradient of the E_BREACH term
			double *	current_energy_E_BREACH		// Current energy of the E_BREACH term
)
{
	// =========================================
	//	Variables

	int i, ind_vertex;

	double * Sigma_vertex = NULL;	// Sum on all disparities using a vertex
	double * Sigma_square_vertex = NULL;	// Sum on all squared disparities using a vertex

	int * ptr_ind_vertex_in_triangle = ind_vertex_in_triangle;
	int * ptr_nb_triangles_using_vertex = NULL;
	double * ptr_D = D, * ptr_grad_D_BREACH = grad_D_BREACH;
	double * ptr_Sigma_vertex = NULL, * ptr_Sigma_square_vertex = NULL;

	// =========================================
	//	Memory Allocation

	Sigma_vertex = (double *)calloc(N_V, sizeof(double));	// Assumed to be allocated with 0.0 values
	Sigma_square_vertex = (double *)calloc(N_V, sizeof(double));


	// =========================================
	//	Compute Sigma_vertex

	for (i = 0; i < 3 * N_T; i++)
	{
		Sigma_vertex[*ptr_ind_vertex_in_triangle] += *ptr_D;
		Sigma_square_vertex[*ptr_ind_vertex_in_triangle] += (*ptr_D) * (*ptr_D);

		ptr_D++;
		ptr_ind_vertex_in_triangle++;
	}



	// =========================================
	//	Compute gradient & current energy

	// Initialization
	*current_energy_E_BREACH = 0.0;

	// Compute gradient
	ptr_D = D;
	ptr_ind_vertex_in_triangle = ind_vertex_in_triangle;

	for (i = 0; i < 3 * N_T; i++)
	{
		*ptr_grad_D_BREACH = 4.0 * (*ptr_D - Sigma_vertex[*ptr_ind_vertex_in_triangle] / nb_triangles_using_vertex[*ptr_ind_vertex_in_triangle]);

		ptr_grad_D_BREACH++;
		ptr_ind_vertex_in_triangle++;
		ptr_D++;
	}

	// Compute energy
	ptr_Sigma_vertex = Sigma_vertex;
	ptr_Sigma_square_vertex = Sigma_square_vertex;
	ptr_nb_triangles_using_vertex = nb_triangles_using_vertex;

	for (i = 0; i < N_V; i++)
	{
		*current_energy_E_BREACH += (*ptr_nb_triangles_using_vertex) * (*ptr_Sigma_square_vertex) - (*ptr_Sigma_vertex) * (*ptr_Sigma_vertex);

		ptr_Sigma_vertex++;
		ptr_Sigma_square_vertex++;
		ptr_nb_triangles_using_vertex++;
	}

	*current_energy_E_BREACH *= 2.0;

	// =========================================
	//	Memory Deallocation

	free(Sigma_vertex);
	free(Sigma_square_vertex);

	return 0;

}










//===========================================================================
///	compute_grad_E_NORMAL
///
/// Compute gradient and current energy of the E_NORMAL term
//===========================================================================
int compute_grad_E_NORMAL(
			int *		ind_triangles_using_edge,	// Indices of the triangles using each edge [size 2*N_E] (max number of triangles using an edge is 2, when less the last indice is -1)
			int			N_E,						// Number of edges
			int			N_T,						// Number of triangles
			double *	A,
			double *	B,
			double *	D,							// Disparity vector
			double *	grad_D_NORMAL,				// Gradient of the E_NORMAL term
			double *	current_energy_E_NORMAL		// Current energy of the E_NORMAL term
)
{
	// =========================================
	//	Variables

	int i, ind_edge, m1, m2;
	double NUM, DENOM, DENOM1, DENOM2;
	
	double A12[9], grad_tmp1[3], grad_tmp2[3];

	int * ptr_ind_triangles_using_edge = ind_triangles_using_edge;
	double * ptr_grad_D_NORMAL = grad_D_NORMAL;
	double * ptr_d1 = NULL, * ptr_d2 = NULL;
	double * ptr_A11 = NULL, *ptr_A22 = NULL;


	// =========================================
	//	Compute gradient & current energy

	// Initialization
	for (i = 0; i < 3 * N_T; i++)
	{
		*ptr_grad_D_NORMAL = 0.0;
		ptr_grad_D_NORMAL++;
	}

	*current_energy_E_NORMAL = 0.0;

	// Main loop
	for (ind_edge = 0; ind_edge < N_E; ind_edge++)
	{
		// Get triangles indices
		m1 = ptr_ind_triangles_using_edge[0];
		m2 = ptr_ind_triangles_using_edge[1];

		// Check if the edge is shared by two triangles (not on image border)
		if (( m1 >= 0 ) && ( m2 >= 0 ))
		{
			ptr_d1 = D + 3 * m1;
			ptr_d2 = D + 3 * m2;

			ptr_A11 = A + 9 * m1;
			ptr_A22 = A + 9 * m2;

			// Compute NUM and DENOM of energy term
			mult_3_3_matrix_by_3_3_transpose_matrix(B + 3 * m1, B + 3 * m2, A12);			// A12 = B1 x B2'
			NUM = mult_1_3_vector_by_3_3_matrix_by_3_1_vector(ptr_d1, A12, D + 3 * m2);		// NUM = d1' x A12 x d2

			// Compute DENOM of energy term
			DENOM1 = mult_1_3_vector_by_3_3_matrix_by_3_1_vector(ptr_d1, ptr_A11, ptr_d1);		// d1' x A11 x d1
			DENOM2 = mult_1_3_vector_by_3_3_matrix_by_3_1_vector(ptr_d2, ptr_A22, ptr_d2);		// d2' x A22 x d2
			DENOM = DENOM1 * DENOM2;


			// Compute gradient for m1
			ptr_grad_D_NORMAL = grad_D_NORMAL + 3 * m1;

			mult_1_3_vector_by_3_3_transpose_matrix(ptr_d2, A12, grad_tmp1);		// grad_tmp1 = d2' x A12'
			mult_1_3_vector_by_3_3_matrix(ptr_d1, ptr_A11, grad_tmp2);			// grad_tmp2 = d1' x A11

			ptr_grad_D_NORMAL[0] -= NUM * (2.0 * DENOM * grad_tmp1[0] - DENOM2 * NUM * grad_tmp2[0]) / (DENOM * DENOM);
			ptr_grad_D_NORMAL[1] -= NUM * (2.0 * DENOM * grad_tmp1[1] - DENOM2 * NUM * grad_tmp2[1]) / (DENOM * DENOM);
			ptr_grad_D_NORMAL[2] -= NUM * (2.0 * DENOM * grad_tmp1[2] - DENOM2 * NUM * grad_tmp2[2]) / (DENOM * DENOM);


			// Compute gradient for m2
			ptr_grad_D_NORMAL = grad_D_NORMAL + 3 * m2;

			mult_1_3_vector_by_3_3_matrix(ptr_d1, A12, grad_tmp1);				// grad_tmp1 = d1' x A12
			mult_1_3_vector_by_3_3_matrix(ptr_d2, ptr_A22, grad_tmp2);			// grad_tmp2 = d2' x A22

			ptr_grad_D_NORMAL[0] -= NUM * (2.0 * DENOM * grad_tmp1[0] - DENOM1 * NUM * grad_tmp2[0]) / (DENOM * DENOM);
			ptr_grad_D_NORMAL[1] -= NUM * (2.0 * DENOM * grad_tmp1[1] - DENOM1 * NUM * grad_tmp2[1]) / (DENOM * DENOM);
			ptr_grad_D_NORMAL[2] -= NUM * (2.0 * DENOM * grad_tmp1[2] - DENOM1 * NUM * grad_tmp2[2]) / (DENOM * DENOM);


			// Compute current energy
			*current_energy_E_NORMAL += 1 - NUM * NUM / DENOM;

		}

		// Update pointers
		ptr_ind_triangles_using_edge += 2;

	}

	return 0;
}
			