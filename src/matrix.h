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

#ifndef matrix3f
typedef GLfloat* matrix3f;
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

BOOL matrix_ortho_set(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far, matrix4f matrix);
BOOL matrix_frustum_set(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat near, GLfloat far, matrix4f matrix);

matrix4f matrix_identity_m4();
matrix3f matrix_identity_m3();
void matrix_identity_m4_set(matrix4f matrix);
void matrix_identity_m3_set(matrix3f matrix);
matrix4f matrix_create_m4(	GLfloat m11, GLfloat m21, GLfloat m31, GLfloat m41,
							GLfloat m12, GLfloat m22, GLfloat m32, GLfloat m42,
							GLfloat m13, GLfloat m23, GLfloat m33, GLfloat m43,
							GLfloat m14, GLfloat m24, GLfloat m34, GLfloat m44);
void matrix_free_m4(const matrix4f matrix);
void matrix_free_m3(const matrix4f matrix);
void matrix_printout_m4(const matrix4f matrix);
void matrix_printout_m3(const matrix3f matrix);

void matrix_get_translation_m4(const matrix4f matrix, GLfloat* x, GLfloat* y, GLfloat* z);
void matrix_get_translation_m3(const matrix3f matrix, GLfloat* x, GLfloat* y);

//count is in number of elements, not components ({x,y,z},{x,y,z} would be 2, not 6)
BOOL matrix_transform_vector(const matrix4f matrix, int dim, int count, const GLfloat* source, GLfloat* dest);
BOOL matrix_transform_normal(const matrix4f matrix, int dim, int count, const GLfloat* source, GLfloat* dest);

matrix4f matrix_scale(GLfloat x, GLfloat y, GLfloat z);
matrix4f matrix_translate(GLfloat x, GLfloat y, GLfloat z);
matrix4f matrix_rotate(GLfloat angle, GLfloat x, GLfloat y, GLfloat z);

void matrix_scale_set(GLfloat x, GLfloat y, GLfloat z, matrix4f matrix);
void matrix_translate_set(GLfloat x, GLfloat y, GLfloat z, matrix4f matrix);
void matrix_rotate_set(GLfloat angle, GLfloat x, GLfloat y, GLfloat z, matrix4f matrix);

GLfloat matrix_determinant(const matrix4f matrix);
matrix4f matrix_transpose_m4(const matrix4f matrix);
matrix3f matrix_transpose_m3(const matrix3f matrix);

BOOL matrix_invert_m4(const matrix4f matrix, matrix4f invert);
BOOL matrix_invert_m3(const matrix3f matrix, matrix3f invert);
BOOL matrix_invert_m4_to_m3(const matrix4f matrix, matrix3f invert);
void matrix_transpose_m4_set(const matrix4f matrix, matrix4f transpose);
void matrix_transpose_m3_set(const matrix3f matrix, matrix3f transpose);

matrix4f matrix_multiply_delete(const matrix4f m1, BOOL freeM1, const matrix4f m2, BOOL freeM2);
matrix4f matrix_multiply(const matrix4f m1, const matrix4f m2);
matrix4f matrix_multiply_disposable(const matrix4f mat, const matrix4f disposable);

#endif
