#ifndef VECTORS_H
#define VECTORS_H



// Returns the norm of the given vector (array).
float vec_norm(float *vec, int size)
{
    int squared_sum = 0;
    for(int i=0; i<size; i++)
        squared_sum += vec[i] * vec[i];
    return sqrt(squared_sum);
}



// Divides each value in the given vector by the given number.
void vec_divide(float *vec, int size, float value)
{
    for(int i=0; i<size; i++)
        vec[i] /= value;
}



// Normalize a vector (make it into a unit vector). If the given vector is
// a zero vector, it is unchanged.
void vec_normalize(float *vec, int size)
{
    float norm = vec_norm(vec, size);
    if(norm > 0)
        vec_divide(vec, size, norm);
}



#endif
