/*
 * matrix.c
 *
 *  Created on: Dec 29, 2012
 *      Author: Vincent Simonetti
 */

#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#include "matrix.h"

matrix4f matrix_ortho(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far)
{
	if(left == right || bottom == top || near == far)
	{
		return NULL;
	}
	GLfloat tx = -(right + left) / (right - left);
	GLfloat ty = -(top + bottom) / (top - bottom);
	GLfloat tz = -(far + near) / (far - near);

	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		ret[0] = 2.0f / (right - left);
		ret[5] = 2.0f / (top - bottom);
		ret[10] = -2.0f / (far - near);

		ret[12] = tx;
		ret[13] = ty;
		ret[14] = tz;

		ret[15] = 1.0f;
	}
	return ret;
}

matrix4f matrix_frustum(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far)
{
	if(left == right || bottom == top || near == far || near < 0.0f || far < 0.0f)
	{
		return NULL;
	}

	GLfloat A = (right + left) / (right - left);
	GLfloat B = (top + bottom) / (top - bottom);
	GLfloat C = -(far + near) / (far - near);
	GLfloat D = -(2.0f * far * near) / (far - near);

	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		ret[0] = (2.0f * near) / (right - left);
		ret[5] = (2.0f * near) / (top - bottom);

		ret[8] = A;
		ret[9] = B;
		ret[10] = C;
		ret[11] = -1.0f;

		ret[14] = D;
	}
	return ret;
}

matrix4f matrix_identiy()
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		ret[0] = 1.0f;
		ret[5] = 1.0f;
		ret[10] = 1.0f;
		ret[15] = 1.0f;
	}
	return ret;
}

matrix4f matrix_create(	GLfloat m11, GLfloat m21, GLfloat m31, GLfloat m41,
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

void matrix_free(const matrix4f matrix)
{
	free(matrix);
}

void matrix_printout(const matrix4f matrix)
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

void matrix_get_translation(const matrix4f matrix, GLfloat* x, GLfloat* y, GLfloat* z)
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

matrix4f matrix_scale(GLfloat x, GLfloat y, GLfloat z)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		ret[0] = x;
		ret[5] = y;
		ret[10] = z;
		ret[15] = 1.0f;
	}
	return ret;
}

matrix4f matrix_translate(GLfloat x, GLfloat y, GLfloat z)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
		ret[0] = 1.0f;
		ret[5] = 1.0f;
		ret[10] = 1.0f;
		ret[15] = 1.0f;

		ret[12] = x;
		ret[13] = y;
		ret[14] = z;
	}
	return ret;
}

matrix4f matrix_rotate(GLfloat angle, GLfloat x, GLfloat y, GLfloat z)
{
	GLfloat* ret = calloc(4 * 4, sizeof(GLfloat));
	if(ret)
	{
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
		GLfloat c = cosf(angle);
		GLfloat s = sinf(angle);

		GLfloat invC = 1.0f - c;

		//Setup the matrix
		ret[0] = ((x * x) * invC) + c;
		ret[1] = (x * y * invC) + (z * s);
		ret[2] = (x * z * invC) - (y * s);

		ret[4] = (x * y * invC) - (z * s);
		ret[5] = ((y * y) * invC) + c;
		ret[6] = (y * z * invC) + (x * s);

		ret[8] = (x * z * invC) + (y * s);
		ret[9] = (y * z * invC) - (x * s);
		ret[10] = ((z * z) * invC) + c;

		ret[15] = 1.0f;
	}
	return ret;
}

void matrix_multiple_internal(GLfloat* m1, GLfloat* m2, GLfloat* product, int sideLen)
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
			matrix_free(m1);
		}
		if(freeM2)
		{
			matrix_free(m2);
		}
	}
	return ret;
}

matrix4f matrix_multiply(const matrix4f m1, const matrix4f m2)
{
	return matrix_multiply_delete(m1, FALSE, m2, FALSE);
}
