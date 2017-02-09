#include <stdio.h>
#include <stdlib.h> 
#include <cmath>

#include "matrix_calculation_functions.h"


// __ Ajouter deux vecteurs

int add_two_vectors(double * A, double * B, int size_vector, double * C)
{
	// C = A + B

	double * ptr_A = A;
	double * ptr_B = B;
	double * ptr_C = C;

	for (int i = 0; i < size_vector; i++)
	{
		*ptr_C = *ptr_A + *ptr_B;

		ptr_A++;
		ptr_B++;
		ptr_C++;
	}

	return 0;
}


// __ Soustraire deux vecteurs

int subtract_two_vectors(double * A, double * B, int size_vector, double * C)
{
	// C = A - B

	double * ptr_A = A;
	double * ptr_B = B;
	double * ptr_C = C;

	for (int i = 0; i < size_vector; i++)
	{
		*ptr_C = *ptr_A - *ptr_B;

		ptr_A++;
		ptr_B++;
		ptr_C++;
	}

	return 0;
}


// __ Copier un vecteur

int copy_vector(double * A, int size_vector, double * B)
{
	// B <- A

	double * ptr_A = A;
	double * ptr_B = B;

	for (int i = 0; i < size_vector; i++)
	{
		*ptr_B = *ptr_A;

		ptr_A++;
		ptr_B++;
	}

	return 0;
}





// __ Multiplication d'un vecteur 1x3 par un vecteur 3x1

double mult_1_3_vector_by_3_1_vector(double * a, double * b)
{
	// c = a x b; 
	double c;

	c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

	return c;
}


// __ Multiplication d'une matrice 3x3 par un vecteur 3x1

int mult_3_1_vector_by_1_3_vector(double * a, double * b, double * C)
{
	// C = a x b; 

	C[0] = a[0] * b[0];
	C[1] = a[0] * b[1];
	C[2] = a[0] * b[2];
	C[3] = a[1] * b[0];
	C[4] = a[1] * b[1];
	C[5] = a[1] * b[2];
	C[6] = a[2] * b[0];
	C[7] = a[2] * b[1];
	C[8] = a[2] * b[2];

	return 0;
}





// __ Multiplication d'une matrice 3x3 par un vecteur 3x1

int mult_3_3_matrix_by_3_1_vector(double * A, double * b, double * c)
{
	// c = A x b;

	c[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
	c[1] = A[3] * b[0] + A[4] * b[1] + A[5] * b[2];
	c[2] = A[6] * b[0] + A[7] * b[1] + A[8] * b[2];

	return 0;
}


// __ Multiplication d'un vecteur 1x3 par une matrice 3x3

int mult_1_3_vector_by_3_3_matrix(double * a, double * B, double * c)
{
	// c = a x B;

	c[0] = a[0] * B[0] + a[1] * B[3] + a[2] * B[6];
	c[1] = a[0] * B[1] + a[1] * B[4] + a[2] * B[7];
	c[2] = a[0] * B[2] + a[1] * B[5] + a[2] * B[8];

	return 0;
}


// __ Multiplication d'un vecteur 1x3 par la transposée d'une matrice 3x3

int mult_1_3_vector_by_3_3_transpose_matrix(double * a, double * B, double * c)
{
	// c = a x B';

	c[0] = a[0] * B[0] + a[1] * B[1] + a[2] * B[2];
	c[1] = a[0] * B[3] + a[1] * B[4] + a[2] * B[5];
	c[2] = a[0] * B[6] + a[1] * B[7] + a[2] * B[8];

	return 0;
}


// __ Multiplication d'un vecteur 1x3 par une matrice 3x3 par un vecteur 3x1
double mult_1_3_vector_by_3_3_matrix_by_3_1_vector(double * a, double * B, double * c)
{
	// d = a x B x c;

	double d;

	d = c[0] * (a[0] * B[0] + a[1] * B[3] + a[2] * B[6]) +
		c[1] * (a[0] * B[1] + a[1] * B[4] + a[2] * B[7]) +
		c[2] * (a[0] * B[2] + a[1] * B[5] + a[2] * B[8]);

	return d;
}





// __ Multiplication de deux matrices 3x3

int mult_two_3_3_matrices(double * A, double * B, double * C)
{
	// C = A x B;

	C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
	C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
	C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];

	C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
	C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
	C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];

	C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
	C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
	C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];

	return 0;
}


// __ Multiplication d'une matrice 3x3 par une matrice 3x3 transposée

int mult_3_3_matrix_by_3_3_transpose_matrix(double * A, double * B, double * C)
{
	// C = A x B';

	C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
	C[1] = A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
	C[2] = A[0] * B[6] + A[1] * B[7] + A[2] * B[8];

	C[3] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
	C[4] = A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
	C[5] = A[3] * B[6] + A[4] * B[7] + A[5] * B[8];

	C[6] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
	C[7] = A[6] * B[3] + A[7] * B[4] + A[8] * B[5];
	C[8] = A[6] * B[6] + A[7] * B[7] + A[8] * B[8];

	return 0;
}


// __ Multiplication d'une matrice 3x3 par sa transposée

int mult_3_3_matrix_by_its_transpose(double * A, double * B)
{
	// B = A x A';

	B[0] = A[0] * A[0] + A[1] * A[1] + A[2] * A[2];
	B[1] = A[0] * A[3] + A[1] * A[4] + A[2] * A[5];
	B[2] = A[0] * A[6] + A[1] * A[7] + A[2] * A[8];

	B[3] = A[3] * A[0] + A[4] * A[1] + A[5] * A[2];
	B[4] = A[3] * A[3] + A[4] * A[4] + A[5] * A[5];
	B[5] = A[3] * A[6] + A[4] * A[7] + A[5] * A[8];

	B[6] = A[6] * A[0] + A[7] * A[1] + A[8] * A[2];
	B[7] = A[6] * A[3] + A[7] * A[4] + A[8] * A[5];
	B[8] = A[6] * A[6] + A[7] * A[7] + A[8] * A[8];

	return 0;
}



// __ Inversion d'une matrice 3x3 quelconque

int invert_3_3_matrix(double * M, double * inv_M)
{
	double det_M;
	double a, b, c, d, e, f, g, h, i;
	double A, B, C, D, E, F, G, H, I;

	// Elements of M
	a = M[0];
	b = M[1];
	c = M[2];
	d = M[3];
	e = M[4];
	f = M[5];
	g = M[6];
	h = M[7];
	i = M[8];

	// Inverse elements
	A = e*i - f*h;
	B = f*g - d*i;
	C = d*h - e*g;
	D = c*h - b*i;
	E = a*i - c*g;
	F = b*g - a*h;
	G = b*f - c*e;
	H = c*d - a*f;
	I = a*e - b*d;

	// Determinant
	det_M = a*A + b*B + c*C;

	// Inversibility test
	if (det_M == 0.0)
	{
		printf("__ Matrice non inversible\n");
		return 1;
	}

	// Inverse matrix
	inv_M[0] = A / det_M;
	inv_M[1] = D / det_M;
	inv_M[2] = G / det_M;
	inv_M[3] = B / det_M;
	inv_M[4] = E / det_M;
	inv_M[5] = H / det_M;
	inv_M[6] = C / det_M;
	inv_M[7] = F / det_M;
	inv_M[8] = I / det_M;

	return 0;
}

