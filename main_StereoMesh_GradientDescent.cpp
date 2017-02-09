#include <stdio.h>
#include <stdlib.h> 
#include <cmath>

#include "matrix_calculation_functions.h"
#include "other_functions_StereoMesh_GradientDescent.h"
#include "energy_StereoMesh_GradientDescent.h"




//===========================================================================
///	StereoMesh_GradientDescent
///
/// StereoMesh using Gradient Descent
//===========================================================================
int StereoMesh_GradientDescent(
			double *	img1,						// Image 1	(RGB)
			double *	img2,						// Image 2	(RGB)
			double *	grad_img2_x,				// Image 2 gradient along x (RGB)		----------- TMP -----------
			int			width_img,					// Image width
			int			height_img,					// Image height
			double *	K,							// Camera calibration matrix of image 1
			double *	B,							// Baseline
			int *		img_label,					// Triangulation labels (segmentation result) [0 to N_T-1]
			int			N_V,						// Number of vertices
			int			N_E,						// Number of edges
			int			N_T,						// Number of triangles
			int *		ind_vertex_in_triangle,		// Indices of the 3 vertices in each triangle [size 3*N_T]
			int *		ind_triangles_using_vertex,	// Indices of the triangles using each vertex [size 6*N_V] (max number of triangles using a vertex is 6, when less the last indices are -1)
			int *		ind_triangles_using_edge,	// Indices of the triangles using each edge [size 2*N_E] (max number of triangles using an edge is 2, when less the last indice is -1)
			double *	D_init,						// Initialization of D
			double *	S,							// S matrices (containing vertices homogeneous coordinates for each triangle)
			double		delta_init,					// Initial step for the gradient descent
			double		lambda_BREACH,				// Weight of E_BREACH
			double		lambda_NORMAL,				// Weight of E_NORMAL
			int			N_ITERmax,					// Max number of iteration for the gradient descent
			double *	final_disparity_map,		// Final disparity_map
			double *	final_img2_interp			// Final interpolation of image 2 from image 1
)
{
	// =========================================
	//	Variables

	int i, N_pixels = height_img * width_img;
	double current_energy_E_DATA, current_energy_E_BREACH, current_energy_E_NORMAL, current_energy;
	double last_energy_E_DATA, last_energy_E_BREACH, last_energy_E_NORMAL, last_energy;

	int * nb_triangles_using_vertex = NULL;

	double * inv_S = NULL, * inv_S_by_pH = NULL;
	double * A = NULL, * B = NULL;
	
	double * D_current = NULL;	// Current disparity vector
	double * D_next = NULL;		// Next disparity vector (if accepted)

	double * grad_D_DATA_current = NULL, * grad_D_BREACH_current = NULL, * grad_D_NORMAL_current = NULL;
	double * grad_D_DATA_next = NULL, * grad_D_BREACH_next = NULL, * grad_D_NORMAL_next = NULL;




	// =========================================
	//	Memory Allocation

	nb_triangles_using_vertex = (int *)calloc(N_V, sizeof(int));
	inv_S = (double *)calloc(9 * N_T, sizeof(double));
	inv_S_by_pH = (double *)calloc(3 * N_pixels, sizeof(double));
	B = (double *)calloc(9 * N_T, sizeof(double));
	A = (double *)calloc(9 * N_T, sizeof(double));

	D_current = (double *)calloc(3 * N_T, sizeof(double));
	D_next = (double *)calloc(3 * N_T, sizeof(double));

	grad_D_DATA_current = (double *)calloc(3 * N_T, sizeof(double));
	grad_D_BREACH_current = (double *)calloc(3 * N_T, sizeof(double));
	grad_D_NORMAL_current = (double *)calloc(3 * N_T, sizeof(double));

	grad_D_DATA_next = (double *)calloc(3 * N_T, sizeof(double));
	grad_D_BREACH_next = (double *)calloc(3 * N_T, sizeof(double));
	grad_D_NORMAL_next = (double *)calloc(3 * N_T, sizeof(double));





	// =========================================
	//	Initialization

	copy_vector(D_init, 3 * N_T, D_current);





	// =========================================
	//	Precomputations

	// __ Inversion des matrices S
	invert_S_matrices(S, N_T, inv_S);

	// __ Calcul des S^(-1) x pH
	compute_inv_S_by_pH(img_label, width_img, height_img, inv_S, N_T, inv_S_by_pH);

	// __ Calcul des matrices B = S^(-1) x K
	compute_B_matrices(inv_S, N_T, K, B);

	// __ Calcul des matrices A = B x B'
	compute_A_matrices(B, N_T, A);

	// __ Calcul de nb_triangles_using_vertex
	compute_nb_triangles_using_vertex(ind_triangles_using_vertex, N_V, nb_triangles_using_vertex);





	// =========================================
	//	First iteration

	// Gradient & energy of the E_DATA term
	compute_grad_E_DATA(img1, img2, grad_img2_x, width_img, height_img, N_pixels, img_label, N_T, inv_S_by_pH, D_current, grad_D_DATA_current, &current_energy_E_DATA);

	// Gradient & energy of the E_BREACH term
	compute_grad_E_BREACH(ind_vertex_in_triangle, ind_triangles_using_vertex, nb_triangles_using_vertex, N_V, N_T, D_current, grad_D_BREACH_current, &current_energy_E_BREACH);

	// Gradient & energy of the E_NORMAL term
	compute_grad_E_NORMAL(ind_triangles_using_edge, N_E, N_T, A, B, D_current, grad_D_NORMAL_current, &current_energy_E_NORMAL);











	// =========================================
	//	Memory Deallocation

	free(nb_triangles_using_vertex);
	free(inv_S);
	free(inv_S_by_pH);
	free(B);
	free(A);

	free(D_current);
	free(D_next);

	free(grad_D_DATA_current);
	free(grad_D_BREACH_current);
	free(grad_D_NORMAL_current);

	free(grad_D_DATA_next);
	free(grad_D_BREACH_next);
	free(grad_D_NORMAL_next);



	return 0;
}