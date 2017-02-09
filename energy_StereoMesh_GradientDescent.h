#pragma once


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
);




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
);


