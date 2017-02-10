#include <stdio.h>
#include <string.h>
#include <vector>

#include "mex.h"

#include "matrix_calculation_functions.cpp"

#include "other_functions_StereoMesh_GradientDescent.cpp"
#include "energy_StereoMesh_GradientDescent.cpp"
#include "main_StereoMesh_GradientDescent.cpp"


//===========================================================================
///	[MEX] StereoMesh_GradientDescent
///
/// StereoMesh using Gradient Descent
//===========================================================================
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	// __ Variables
	int height_img, width_img;
    
	int * dims_img1 = NULL;
	int * dims_img2 = NULL;
	int * dims_grad_img2_x = NULL;

	int dims_img_out[2];

	// __ Entrées
	int N_V;		// Nombre de sommets
	int N_E;		// Nombre d'arêtes
	int N_T;		// Nombre de triangles
	int	N_ITERmax;	// Max number of iteration for the gradient descent

	double Baseline;				// Baseline
	double delta_init;				// Initial step for the gradient descent
	double lambda_BREACH;			// Weight of E_BREACH
	double lambda_NORMAL;			// Weight of E_NORMAL

	int * img_label = NULL;
	
	int * ind_vertex_in_triangle = NULL;		// Indices of the 3 vertices in each triangle [size 3*N_T]
	int * ind_triangles_using_vertex = NULL;	// Indices of the triangles using each vertex [size 6*N_V] (max number of triangles using a vertex is 6, when less the last indices are -1)
	int * ind_triangles_using_edge = NULL;		// Indices of the triangles using each edge [size 2*N_E] (max number of triangles using an edge is 2, when less the last indice is -1)

	double * S = NULL;
	double * D_init = NULL;
	
	double * img1 = NULL, * img2 = NULL;
	double * grad_img2_x = NULL;

	double * K = NULL;

	// __ Sorties
	double * final_disparity_map = NULL;		// Final disparity_map
	double * final_img2_interp = NULL;


    // __ Récupération des entrées
    img1 = (double *)mxGetPr(prhs[0]);
	img2 = (double *)mxGetPr(prhs[1]);
	grad_img2_x = (double *)mxGetPr(prhs[2]);
	K = (double *)mxGetPr(prhs[3]);
	Baseline = (double)mxGetScalar(prhs[4]);
	img_label = (int *)mxGetPr(prhs[5]);
	N_V = (int)mxGetScalar(prhs[6]);
	N_E = (int)mxGetScalar(prhs[7]);
	N_T = (int)mxGetScalar(prhs[8]);
	ind_vertex_in_triangle = (int *)mxGetPr(prhs[9]);
	ind_triangles_using_vertex = (int *)mxGetPr(prhs[10]);
	ind_triangles_using_edge = (int *)mxGetPr(prhs[11]);
	D_init = (double *)mxGetPr(prhs[12]);
	S = (double *)mxGetPr(prhs[13]);
	delta_init = (double)mxGetScalar(prhs[14]);
	lambda_BREACH = (double)mxGetScalar(prhs[15]);
	lambda_NORMAL = (double)mxGetScalar(prhs[16]);
	N_ITERmax = (int)mxGetScalar(prhs[17]);

    // __ Dimensions des entrées
	dims_img1 = (int *)mxGetDimensions(prhs[0]);
	dims_img2 = (int *)mxGetDimensions(prhs[1]);
	dims_grad_img2_x = (int *)mxGetDimensions(prhs[2]);
	
	// __ Tests sur les entrées
	if ((mxGetNumberOfDimensions(prhs[0]) != 3) || (mxGetNumberOfDimensions(prhs[1]) != 3) || (mxGetNumberOfDimensions(prhs[2]) != 3))
	{
		mexPrintf("mexERROR - Une des images n'est pas de dimension 3..................\n");
		return;
	}

	if ((dims_img1[0] != dims_img2[0]) || (dims_img1[0] != dims_grad_img2_x[0]) || (dims_img1[1] != dims_img2[1]) || (dims_img1[1] != dims_grad_img2_x[1]))
	{
		mexPrintf("mexERROR - Dimensions des images différentes..................\n");
		return;
	}
	
	if ((dims_img1[2] != 3) || (dims_img2[2] != 3) || (dims_grad_img2_x[2] != 3))
	{
		mexPrintf("mexERROR - Images RGB attendues..................\n");
		return;
	}
	
    width_img = dims_img1[0];
    height_img = dims_img1[1];

//	mexPrintf("width_img : %d\n", width_img);
//	mexPrintf("height_img : %d\n", height_img);

	// __ Création des sorties
	dims_img_out[0] = width_img;
	dims_img_out[1] = height_img;

	plhs[0] = mxCreateNumericArray(2, dims_img_out, mxDOUBLE_CLASS, mxREAL);
	final_disparity_map = (double *)mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericArray(2, dims_img_out, mxDOUBLE_CLASS, mxREAL);
	final_img2_interp = (double *)mxGetPr(plhs[1]);

	

    // __ Programme
	StereoMesh_GradientDescent(img1, img2, grad_img2_x, width_img, height_img, K, Baseline, img_label, N_V, N_E, N_T, ind_vertex_in_triangle, ind_triangles_using_vertex,
									ind_triangles_using_edge, D_init, S, delta_init, lambda_BREACH, lambda_NORMAL, N_ITERmax, final_disparity_map, final_img2_interp);

}
