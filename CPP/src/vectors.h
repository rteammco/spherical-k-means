/* File: vectors.h
 *
 * This file contains basic vector manipulation functions used by the
 * spherical K-means algorithm.
 */

#ifndef VECTORS_H
#define VECTORS_H


// Returns the norm of the given vector (array).
float vec_norm(float *vec, int size);

// Returns the sum of each component of the the given vector.
float vec_sum(float *vec, int size);

// Returns the dot product of the two given vectors.
float vec_dot(float *vec1, float *vec2, int size);

// [return] Returns a new vector that is the sum of the given list of vectors.
float* vec_sum(float **vecs, int size, int num_vecs);

// [return] Scales each value in the given vector by the given power,
// and returns the resulting vector.
float* vec_pow_new(float *vec, int size, float power);

// [return] Returns a new zero-vector of the given size.
float* vec_zeros(int size);

// [in-place] Adds the second vector to the first one.
void vec_add(float *vec1, float *vec2, int size);

// [in-place] Multiplies each value in the given vector by the given number.
void vec_multiply(float *vec, int size, float value);

// [in-place] Divides each value in the given vector by the given number.
void vec_divide(float *vec, int size, float value);

// [in-place] Scales each value in the given vector by the given power.
void vec_pow(float *vec, int size, float power);

// [in-place] Normalize a vector (make it into a unit vector). If the given
// vector is a zero vector, it is unchanged.
void vec_normalize(float *vec, int size);


#endif
