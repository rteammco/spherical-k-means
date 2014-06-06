#ifndef VECTORS_H
#define VECTORS_H

#include <math.h>



// Returns the norm of the given vector (array).
float vec_norm(float *vec, int size)
{
    float squared_sum = 0;
    for(int i=0; i<size; i++)
        squared_sum += vec[i] * vec[i];
    return sqrt(squared_sum);
}



// Returns the sum of each component of the the given vector.
float vec_sum(float *vec, int size)
{
    float sum = 0;
    for(int i=0; i<size; i++)
        sum += vec[i];
    return sum;
}



// Returns the dot product of the two given vectors.
float vec_dot(float *vec1, float *vec2, int size)
{
    float dotp = 0;
    for(int i=0; i<size; i++)
        dotp += vec1[i] * vec2[i];
    return dotp;
}



// [return] Returns a new vector that is the sum of the given list of vectors.
float* vec_sum(float **vecs, int size, int num_vecs)
{
    float *sum = new float[size];
    for(int i=0; i<size; i++) {
        sum[i] = 0;
        for(int j=0; j<num_vecs; j++)
            sum[i] += vecs[j][i];
    }
    return sum;
}



// [return] Scales each value in the given vector by the given power,
// and returns the resulting vector.
float* vec_pow_new(float *vec, int size, float power)
{
    float *new_vec = new float[size];
    for(int i=0; i<size; i++) {
        if(vec[i] > 0)
            new_vec[i] = pow(vec[i], power);
        else
            new_vec[i] = 0;
    }
    return new_vec;
}



// [in-place] Multiplies each value in the given vector by the given number.
void vec_multiply(float *vec, int size, float value)
{
    for(int i=0; i<size; i++)
        vec[i] *= value;
}



// [in-place] Divides each value in the given vector by the given number.
void vec_divide(float *vec, int size, float value)
{
    for(int i=0; i<size; i++)
        vec[i] /= value;
}



// [in-place] Scales each value in the given vector by the given power.
void vec_pow(float *vec, int size, float power)
{
    for(int i=0; i<size; i++)
        vec[i] = pow(vec[i], power);
}



// [in-place] Normalize a vector (make it into a unit vector). If the given
// vector is a zero vector, it is unchanged.
void vec_normalize(float *vec, int size)
{
    float norm = vec_norm(vec, size);
    if(norm > 0)
        vec_divide(vec, size, norm);
}



#endif
