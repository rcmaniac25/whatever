/*
 * matrix.c
 *
 *  Created on: Dec 29, 2012
 *      Author: Vincent Simonetti
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <math.h>

#include "matrix.h"
#include "vvector.h"

//Nothing in here is optimized, so it would be easier to understand.

matrix4f matrix_ortho(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		if(!matrix_ortho_set(left, right, bottom, top, near, far, ret))
		{
			free(ret);
			return NULL;
		}
	}
	return ret;
}

matrix4f matrix_frustum(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		if(!matrix_frustum_set(left, right, bottom, top, near, far, ret))
		{
			free(ret);
			return NULL;
		}
	}
	return ret;
}

BOOL matrix_ortho_set(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far, matrix4f matrix)
{
	if(!matrix || left == right || bottom == top || near == far)
	{
		return FALSE;
	}
	GLfloat tx = -(right + left) / (right - left);
	GLfloat ty = -(top + bottom) / (top - bottom);
	GLfloat tz = -(far + near) / (far - near);

	memset(matrix, 0, sizeof(GLfloat) * 4 * 4);

	matrix[0] = 2.0f / (right - left);
	matrix[5] = 2.0f / (top - bottom);
	matrix[10] = -2.0f / (far - near);

	matrix[12] = tx;
	matrix[13] = ty;
	matrix[14] = tz;

	matrix[15] = 1.0f;

	return TRUE;
}

BOOL matrix_frustum_set(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far, matrix4f matrix)
{
	if(!matrix || left == right || bottom == top || near == far || near < 0.0f || far < 0.0f)
	{
		return FALSE;
	}

	GLfloat A = (right + left) / (right - left);
	GLfloat B = (top + bottom) / (top - bottom);
	GLfloat C = -(far + near) / (far - near);
	GLfloat D = -(2.0f * far * near) / (far - near);

	memset(matrix, 0, sizeof(GLfloat) * 4 * 4);

	matrix[0] = (2.0f * near) / (right - left);
	matrix[5] = (2.0f * near) / (top - bottom);

	matrix[8] = A;
	matrix[9] = B;
	matrix[10] = C;
	matrix[11] = -1.0f;

	matrix[14] = D;

	return TRUE;
}

matrix4f matrix_identity_m4()
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	matrix_identity_m4_set(ret);
	return ret;
}

matrix3f matrix_identity_m3()
{
	GLfloat* ret = calloc(3 * 3, sizeof(GLfloat));
	matrix_identity_m3_set(ret);
	return ret;
}

void matrix_identity_m4_set(matrix4f matrix)
{
	if(matrix)
	{
		memset(matrix, 0, sizeof(GLfloat) * 4 * 4);

		matrix[0] = 1.0f;
		matrix[5] = 1.0f;
		matrix[10] = 1.0f;
		matrix[15] = 1.0f;
	}
}

void matrix_identity_m3_set(matrix3f matrix)
{
	if(matrix)
	{
		memset(matrix, 0, sizeof(GLfloat) * 3 * 3);

		matrix[0] = 1.0f;
		matrix[5] = 1.0f;
		matrix[10] = 1.0f;
	}
}

matrix4f matrix_create_m4(	GLfloat m11, GLfloat m21, GLfloat m31, GLfloat m41,
							GLfloat m12, GLfloat m22, GLfloat m32, GLfloat m42,
							GLfloat m13, GLfloat m23, GLfloat m33, GLfloat m43,
							GLfloat m14, GLfloat m24, GLfloat m34, GLfloat m44)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		ret[0] = m11;
		ret[1] = m12;
		ret[2] = m13;
		ret[3] = m14;

		ret[4] = m21;
		ret[5] = m22;
		ret[6] = m23;
		ret[7] = m24;

		ret[8] = m31;
		ret[9] = m32;
		ret[10] = m33;
		ret[11] = m34;

		ret[12] = m41;
		ret[13] = m42;
		ret[14] = m43;
		ret[15] = m44;
	}
	return ret;
}

void matrix_free_m4(const matrix4f matrix)
{
	free(matrix);
}

void matrix_free_m3(const matrix3f matrix)
{
	free(matrix);
}

void matrix_printout_m4(const matrix4f matrix)
{
	if(matrix)
	{
		fprintf(stderr, "%f, %f, %f, %f\n", matrix[0], matrix[4], matrix[8], matrix[12]);
		fprintf(stderr, "%f, %f, %f, %f\n", matrix[1], matrix[5], matrix[9], matrix[13]);
		fprintf(stderr, "%f, %f, %f, %f\n", matrix[2], matrix[6], matrix[10], matrix[14]);
		fprintf(stderr, "%f, %f, %f, %f\n", matrix[3], matrix[7], matrix[11], matrix[15]);
		fprintf(stderr, "\n");
	}
}

void matrix_printout_m3(const matrix3f matrix)
{
	if(matrix)
	{
		fprintf(stderr, "%f, %f, %f\n", matrix[0], matrix[3], matrix[6]);
		fprintf(stderr, "%f, %f, %f\n", matrix[1], matrix[4], matrix[7]);
		fprintf(stderr, "%f, %f, %f\n", matrix[2], matrix[5], matrix[8]);
		fprintf(stderr, "\n");
	}
}

void matrix_get_translation_m4(const matrix4f matrix, GLfloat* x, GLfloat* y, GLfloat* z)
{
	if(matrix)
	{
		if(x)
		{
			*x = matrix[12];
		}
		if(y)
		{
			*y = matrix[13];
		}
		if(z)
		{
			*z = matrix[14];
		}
	}
}

void matrix_get_translation_m3(const matrix3f matrix, GLfloat* x, GLfloat* y)
{
	if(matrix)
	{
		if(x)
		{
			*x = matrix[6];
		}
		if(y)
		{
			*y = matrix[7];
		}
	}
}

void matrix_transform_internal(const GLfloat* matrix, int sideLen, int dim, int count, BOOL translate, const GLfloat* source, GLfloat* dest)
{
	int i;
	int c;
	int m;
	for(i = 0; i < count; i++)
	{
		for(c = 0; c < dim; c++)
		{
			int destPos = i * dim + c;
			dest[destPos] = 0.0f;
			for(m = 0; m < dim; m++)
			{
				dest[destPos] += source[i * dim + c] * matrix[m * dim + c];
			}
			if(translate && dim < sideLen)
			{
				dest[destPos] += matrix[(sideLen - 1) * sideLen + c];
			}
		}
	}
}

BOOL matrix_transform_vector(const matrix4f matrix, int dim, int count, const GLfloat* source, GLfloat* dest)
{
	if(!matrix || dim < 1 || dim > 4 || !source || !dest)
	{
		return FALSE;
	}
	matrix_transform_internal(matrix, 4, dim, count, TRUE, source, dest);
	return TRUE;
}

BOOL matrix_transform_normal(const matrix4f matrix, int dim, int count, const GLfloat* source, GLfloat* dest)
{
	if(!matrix || dim < 1 || dim > 4 || !source || !dest)
	{
		return FALSE;
	}
	matrix_transform_internal(matrix, 4, dim, count, FALSE, source, dest);
	return TRUE;
}

matrix4f matrix_scale(GLfloat x, GLfloat y, GLfloat z)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	matrix_scale_set(x, y, z, ret);
	return ret;
}

matrix4f matrix_translate(GLfloat x, GLfloat y, GLfloat z)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	matrix_translate_set(x, y, z, ret);
	return ret;
}

matrix4f matrix_rotate(GLfloat angle, GLfloat x, GLfloat y, GLfloat z)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	matrix_rotate_set(angle, x, y, z, ret);
	return ret;
}

void matrix_scale_set(GLfloat x, GLfloat y, GLfloat z, matrix4f matrix)
{
	if(matrix)
	{
		memset(matrix, 0, sizeof(GLfloat) * 4 * 4);

		matrix[0] = x;
		matrix[5] = y;
		matrix[10] = z;
		matrix[15] = 1.0f;
	}
}

void matrix_translate_set(GLfloat x, GLfloat y, GLfloat z, matrix4f matrix)
{
	if(matrix)
	{
		memset(matrix, 0, sizeof(GLfloat) * 4 * 4);

		matrix[0] = 1.0f;
		matrix[5] = 1.0f;
		matrix[10] = 1.0f;
		matrix[15] = 1.0f;

		matrix[12] = x;
		matrix[13] = y;
		matrix[14] = z;
	}
}

void matrix_rotate_set(GLfloat angle, GLfloat x, GLfloat y, GLfloat z, matrix4f matrix)
{
	if(matrix)
	{
		//Convert the angle to radians
		GLfloat rad = angle * (M_PI / 180.0f);

		//Normalize input
		GLfloat len = sqrtf((x * x) + (y * y) + (z * z));
		if(len != 0.0f)
		{
			len = 1.0f / len;
			x *= len;
			y *= len;
			z *= len;
		}

		//Calculate sin/cos
		GLfloat c = cosf(rad);
		GLfloat s = sinf(rad);

		GLfloat invC = 1.0f - c;

		//Setup the matrix
		matrix[0] = ((x * x) * invC) + c;
		matrix[1] = (x * y * invC) + (z * s);
		matrix[2] = (x * z * invC) - (y * s);
		matrix[3] = 0.0f;

		matrix[4] = (x * y * invC) - (z * s);
		matrix[5] = ((y * y) * invC) + c;
		matrix[6] = (y * z * invC) + (x * s);
		matrix[7] = 0.0f;

		matrix[8] = (x * z * invC) + (y * s);
		matrix[9] = (y * z * invC) - (x * s);
		matrix[10] = ((z * z) * invC) + c;
		matrix[11] = 0.0f;

		matrix[12] = 0.0f;
		matrix[13] = 0.0f;
		matrix[14] = 0.0f;
		matrix[15] = 1.0f;
	}
}

#define MATRIX4F_TO_MAT4X4(mat,m44) \
	m44[0][0] = mat[0 * 4 + 0]; \
	m44[0][1] = mat[0 * 4 + 1]; \
	m44[0][2] = mat[0 * 4 + 2]; \
	m44[0][3] = mat[0 * 4 + 3]; \
	m44[1][0] = mat[1 * 4 + 0]; \
	m44[1][1] = mat[1 * 4 + 1]; \
	m44[1][2] = mat[1 * 4 + 2]; \
	m44[1][3] = mat[1 * 4 + 3]; \
	m44[2][0] = mat[2 * 4 + 0]; \
	m44[2][1] = mat[2 * 4 + 1]; \
	m44[2][2] = mat[2 * 4 + 2]; \
	m44[2][3] = mat[2 * 4 + 3]; \
	m44[3][0] = mat[3 * 4 + 0]; \
	m44[3][1] = mat[3 * 4 + 1]; \
	m44[3][2] = mat[3 * 4 + 2]; \
	m44[3][3] = mat[3 * 4 + 3];

#define MATRIX3F_TO_MAT3X3(mat,m33) \
	m33[0][0] = mat[0 * 3 + 0]; \
	m33[0][1] = mat[0 * 3 + 1]; \
	m33[0][2] = mat[0 * 3 + 2]; \
	m33[1][0] = mat[1 * 3 + 0]; \
	m33[1][1] = mat[1 * 3 + 1]; \
	m33[1][2] = mat[1 * 3 + 2]; \
	m33[2][0] = mat[2 * 3 + 0]; \
	m33[2][1] = mat[2 * 3 + 1]; \
	m33[2][2] = mat[2 * 3 + 2];

#define MAT4X4_TO_MATRIX4F(m44,mat) \
	mat[0 * 4 + 0] = m44[0][0]; \
	mat[0 * 4 + 1] = m44[0][1]; \
	mat[0 * 4 + 2] = m44[0][2]; \
	mat[0 * 4 + 3] = m44[0][3]; \
	mat[1 * 4 + 0] = m44[1][0]; \
	mat[1 * 4 + 1] = m44[1][1]; \
	mat[1 * 4 + 2] = m44[1][2]; \
	mat[1 * 4 + 3] = m44[1][3]; \
	mat[2 * 4 + 0] = m44[2][0]; \
	mat[2 * 4 + 1] = m44[2][1]; \
	mat[2 * 4 + 2] = m44[2][2]; \
	mat[2 * 4 + 3] = m44[2][3]; \
	mat[3 * 4 + 0] = m44[3][0]; \
	mat[3 * 4 + 1] = m44[3][1]; \
	mat[3 * 4 + 2] = m44[3][2]; \
	mat[3 * 4 + 3] = m44[3][3];

#define MAT3X3_TO_MATRIX3F(m33,mat) \
	mat[0 * 3 + 0] = m33[0][0]; \
	mat[0 * 3 + 1] = m33[0][1]; \
	mat[0 * 3 + 2] = m33[0][2]; \
	mat[1 * 3 + 0] = m33[1][0]; \
	mat[1 * 3 + 1] = m33[1][1]; \
	mat[1 * 3 + 2] = m33[1][2]; \
	mat[2 * 3 + 0] = m33[2][0]; \
	mat[2 * 3 + 1] = m33[2][1]; \
	mat[2 * 3 + 2] = m33[2][2];

#define MAT4X4_TO_MATRIX3F(m44,mat) \
	mat[0 * 3 + 0] = m44[0][0]; \
	mat[0 * 3 + 1] = m44[0][1]; \
	mat[0 * 3 + 2] = m44[0][2]; \
	mat[1 * 3 + 0] = m44[1][0]; \
	mat[1 * 3 + 1] = m44[1][1]; \
	mat[1 * 3 + 2] = m44[1][2]; \
	mat[2 * 3 + 0] = m44[2][0]; \
	mat[2 * 3 + 1] = m44[2][1]; \
	mat[2 * 3 + 2] = m44[2][2];

GLfloat matrix_determinant(const matrix4f matrix)
{
	if(matrix)
	{
		GLfloat mat[4][4];
		GLfloat d = 0.0;

		MATRIX4F_TO_MAT4X4(matrix,mat)
		DETERMINANT_4X4(d,mat)

		return d;
	}
	return INFINITY;
}

matrix4f matrix_transpose_m4(const matrix4f matrix)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	matrix_transpose_m4_set(matrix, ret);
	return ret;
}

matrix3f matrix_transpose_m3(const matrix3f matrix)
{
	GLfloat* ret = calloc(3 * 3, sizeof(GLfloat));
	matrix_transpose_m3_set(matrix, ret);
	return ret;
}

BOOL matrix_invert_m4(const matrix4f matrix, matrix4f invert)
{
	if(matrix && invert && matrix != invert)
	{
		GLfloat in[4][4];
		GLfloat out[4][4];
		GLfloat d = 0.0;

		MATRIX4F_TO_MAT4X4(matrix,in);
		INVERT_4X4(out,d,in);
		if(d >= 1e-6)
		{
			MAT4X4_TO_MATRIX4F(out,invert);
			return TRUE;
		}
	}
	return FALSE;
}

BOOL matrix_invert_m3(const matrix3f matrix, matrix3f invert)
{
	if(matrix && invert && matrix != invert)
	{
		GLfloat in[3][3];
		GLfloat out[3][3];
		GLfloat d = 0.0;

		MATRIX3F_TO_MAT3X3(matrix,in);
		INVERT_3X3(out,d,in);
		if(d >= 1e-6)
		{
			MAT3X3_TO_MATRIX3F(out,invert);
			return TRUE;
		}
	}
	return FALSE;
}

BOOL matrix_invert_m4_to_m3(const matrix4f matrix, matrix3f invert)
{
	if(matrix && invert && matrix != invert)
	{
		GLfloat in[4][4];
		GLfloat out[4][4];
		GLfloat d = 0.0;

		MATRIX4F_TO_MAT4X4(matrix,in);
		INVERT_4X4(out,d,in);
		if(d >= 1e-6)
		{
			MAT4X4_TO_MATRIX3F(out,invert);
			return TRUE;
		}
	}
	return FALSE;
}

void matrix_transpose_internal(const GLfloat* matrix, GLfloat* transpose, int sideLen)
{
	int i, j;
	for(i = 0; i < sideLen; i++)
	{
		for(j = 0; j < sideLen; j++)
		{
			transpose[j * sideLen + i] = matrix[i * sideLen + j];
		}
	}
}

void matrix_transpose_m4_set(const matrix4f matrix, matrix4f transpose)
{
	if(matrix && transpose && matrix != transpose)
	{
		matrix_transpose_internal(matrix, transpose, 4);
	}
}

void matrix_transpose_m3_set(const matrix3f matrix, matrix3f transpose)
{
	if(matrix && transpose && matrix != transpose)
	{
		matrix_transpose_internal(matrix, transpose, 3);
	}
}

void matrix_multiple_internal(const GLfloat* m1, const GLfloat* m2, GLfloat* product, int sideLen)
{
	int i, j, k;

	//Not efficient but simple
	for(i = 0; i < sideLen; i++)
	{
		for(j = 0; j < sideLen; j++)
		{
			for(k = 0; k < sideLen; k++)
			{
				product[i * sideLen + j] += m1[i * sideLen + k] * m2[k * sideLen + j];
			}
		}
	}
}

//m1 * m2
matrix4f matrix_multiply_delete(const matrix4f m1, BOOL freeM1, const matrix4f m2, BOOL freeM2)
{
	if(!m1 || !m2)
	{
		return NULL;
	}
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		//Not efficient but simple
		matrix_multiple_internal(m2, m1, ret, 4);

		if(freeM1)
		{
			matrix_free_m4(m1);
		}
		if(freeM2)
		{
			matrix_free_m4(m2);
		}
	}
	return ret;
}

matrix4f matrix_multiply(const matrix4f m1, const matrix4f m2)
{
	return matrix_multiply_delete(m1, FALSE, m2, FALSE);
}

matrix4f matrix_multiply_disposable(const matrix4f mat, const matrix4f disposable)
{
	return matrix_multiply_delete(mat, FALSE, disposable, TRUE);
}
