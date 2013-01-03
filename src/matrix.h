/*
 * matrix.h
 *
 *  Created on: Dec 29, 2012
 *      Author: Vincent Simonetti (Rcmaniac25)
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#ifndef GLfloat
typedef float GLfloat;
#endif

#ifndef matrix4f
typedef GLfloat* matrix4f;
#endif

#ifndef BOOL
#define BOOL int
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

matrix4f matrix_ortho(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far);
matrix4f matrix_frustum(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far);

matrix4f matrix_identiy();
matrix4f matrix_create(	GLfloat m11, GLfloat m21, GLfloat m31, GLfloat m41,
						GLfloat m12, GLfloat m22, GLfloat m32, GLfloat m42,
						GLfloat m13, GLfloat m23, GLfloat m33, GLfloat m43,
						GLfloat m14, GLfloat m24, GLfloat m34, GLfloat m44);
void matrix_free(const matrix4f matrix);
void matrix_printout(const matrix4f matrix);

void matrix_get_translation(const matrix4f matrix, GLfloat* x, GLfloat* y, GLfloat* z);

//count is in number of elements, not components ({x,y,z},{x,y,z} would be 2, not 6)
BOOL matrix_transform_vector(const matrix4f matrix, int dim, int count, const GLfloat* source, GLfloat* dest);
BOOL matrix_transform_normal(const matrix4f matrix, int dim, int count, const GLfloat* source, GLfloat* dest);

matrix4f matrix_scale(GLfloat x, GLfloat y, GLfloat z);
matrix4f matrix_translate(GLfloat x, GLfloat y, GLfloat z);
matrix4f matrix_rotate(GLfloat angle, GLfloat x, GLfloat y, GLfloat z);

matrix4f matrix_multiply_delete(const matrix4f m1, BOOL freeM1, const matrix4f m2, BOOL freeM2);
matrix4f matrix_multiply(const matrix4f m1, const matrix4f m2);

#endif
