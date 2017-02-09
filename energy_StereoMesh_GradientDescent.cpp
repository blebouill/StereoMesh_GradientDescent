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
				*current_energy_E_DATA += delta_img21_R*delta_img21_R + delta_img21_G*delta_img21_G + delta21_grad_B*delta21_grad_B;

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
			ptr_grad_D_DATA[0] /= (double)(*ptr_norm_grad_nb_pix);
			ptr_grad_D_DATA[1] /= (double)(*ptr_norm_grad_nb_pix);
			ptr_grad_D_DATA[2] /= (double)(*ptr_norm_grad_nb_pix);
		}
		ptr_norm_grad_nb_pix++;
		ptr_grad_D_DATA += 3;
	}



	// =========================================
	//	Memory Deallocation

	free(norm_grad_nb_pix);

	return 0;
}