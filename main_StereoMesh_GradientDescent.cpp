#include <stdio.h>
#include <stdlib.h> 
#include <cmath>

#include "matrix_calculation_functions.h"
#include "other_functions_StereoMesh_GradientDescent.h"




//===========================================================================
///	StereoMesh_GradientDescent
///
/// StereoMesh using Gradient Descent
//===========================================================================
int StereoMesh_GradientDescent(
			double *	img1,					// Image 1	(RGB)
			double *	img2,					// Image 2	(RGB)
			double *	grad_img2_x,			// Image 2 gradient along x (RGB)		----------- TMP -----------
			int			width_img,				// Image width
			int			height_img,				// Image height
			double *	K,						// Camera calibration matrix of image 1
			double *	B,						// Baseline
			int *		img_label,				// Triangulation labels (segmentation result) [0 to N_T-1]
			int			N_V,					// Number of vertices
			int			N_E,					// Number of edges
			int			N_T,					// Number of triangles
			double *	S,						// S matrices (containing vertices homogeneous coordinates for each triangle)
			int			N_ITERmax,				// Max number of iteration for the gradient descent
			double *	final_disparity_map,	// Final disparity_map
			double *	final_img2_interp		// Final interpolation of image 2 from image 1
)
{
	// =========================================
	//	Variables

	int N_pixels = height_img * width_img;

	double * inv_S = NULL, * inv_S_by_pH = NULL;
	double * A = NULL, * B = NULL;
	double * D = NULL;	// Disparity vector






	// =========================================
	//	Memory Allocation

	inv_S = (double *)calloc(9 * N_T, sizeof(double));
	inv_S_by_pH = (double *)calloc(3 * N_pixels, sizeof(double));
	B = (double *)calloc(9 * N_T, sizeof(double));
	A = (double *)calloc(9 * N_T, sizeof(double));
	D = (double *)calloc(3 * N_T, sizeof(double));





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








	// =========================================
	//	Memory Deallocation

	free(inv_S);
	free(inv_S_by_pH);
	free(B);
	free(A);
	free(D);

	return 0;
}