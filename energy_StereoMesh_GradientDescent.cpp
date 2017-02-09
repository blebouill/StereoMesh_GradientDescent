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
			double *	grad_img2_x,			// Image 2 gradient along x (RGB)		----------- TMP -----------
			int			width_img,				// Image width
			int			height_img,				// Image height
			int			N_pixels,				// N_pixels = height_img x width_img;
			int *		img_label,				// Triangulation labels (segmentation result) [0 to N_T-1]
			int			N_T,					// Number of triangles
			double *	inv_S_by_pH,			// S^(-1) x pH
			double *	D,						// Disparity vector
			double *	grad_D_DATA				// Gradient of the E_DATA term
)
{
	// =========================================
	//	Variables

	int i, x, y, m;
	double x_img2;

	int * norm_grad_nb_pix = NULL;

	int * ptr_img_label = img_label;

	double * ptr_inv_S_by_pH = inv_S_by_pH;

	double * ptr_img1_R = img1;
	double * ptr_img1_G = img1 + N_pixels;
	double * ptr_img1_B = img1 + 2 * N_pixels;

	double * ptr_img2_R = img1;
	double * ptr_img2_G = img1 + N_pixels;
	double * ptr_img2_B = img1 + 2 * N_pixels;

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
		grad_D_DATA[i] = 0.0;

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

			}





			// Update pointers
			ptr_img_label++;

			ptr_inv_S_by_pH++;

			ptr_img1_R++;
			ptr_img1_G++;
			ptr_img1_B++;

			ptr_img2_R++;
			ptr_img2_G++;
			ptr_img2_B++;

			ptr_grad_img2_x_R++;
			ptr_grad_img2_x_G++;
			ptr_grad_img2_x_B++;
		}
	}



	// =========================================
	//	Memory Deallocation

	free(norm_grad_nb_pix);

	return 0;
}