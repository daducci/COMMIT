#ifndef __OPENGL_UTILS_H__
#define __OPENGL_UTILS_H__

#include <algorithm>

#include "VECTOR.h"
typedef VECTOR<GLint>		Vec3Di;
typedef VECTOR<GLfloat>		Vec3Df;


namespace OPENGL_utils
{

void identity(GLfloat* result)
{
    for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
        if (i==j) result[4*i+j]=1; else result[4*i+j]=0;
}


void matXMat(GLfloat* m, GLfloat* m1, GLfloat* result)
{
    for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
        result[4*i+j]=0;
        for (int t=0; t<4; t++)
            result[4*i+j]=result[4*i+j]+m[4*i+t]*m1[4*t+j];
    }
}


void rotateZ(GLfloat* m, GLfloat ang, GLfloat* result)
{
    static GLfloat matrix[16];

    for (int i=0; i<16 ; i++) matrix[i] = 0;
    matrix[0]  = cos(ang/180*3.1415);
    matrix[5]  = cos(ang/180*3.1415);
    matrix[1]  = -sin(ang/180*3.1415);
    matrix[4]  = sin(ang/180*3.1415);
    matrix[10] = 1;
    matrix[15] = 1;
    matXMat(matrix,m,result);
}


void rotateY(GLfloat* m, GLfloat ang, GLfloat* result)
{
    static GLfloat matrix[16];

    for (int i=0; i<16 ; i++) matrix[i] = 0;
    matrix[0]  = cos(ang/180*3.1415);
    matrix[10] = cos(ang/180*3.1415);
    matrix[8]  = -sin(ang/180*3.1415);
    matrix[2]  = sin(ang/180*3.1415);
    matrix[5]  = 1;
    matrix[15] = 1;
    matXMat(matrix,m,result);
}


void rotateX(GLfloat* m, GLfloat ang, GLfloat* result)
{
    static GLfloat matrix[16];

    for (int i=0; i<16 ; i++) matrix[i] = 0;
    matrix[5]  = cos(ang/180*3.1415);
    matrix[10] = cos(ang/180*3.1415);
    matrix[6]  = -sin(ang/180*3.1415);
    matrix[9]  = sin(ang/180*3.1415);
    matrix[0]  = 1;
    matrix[15] = 1;
    matXMat(matrix,m,result);
}


void translate(GLfloat* m, GLfloat x,GLfloat y,GLfloat z, GLfloat* result)
{
    static GLfloat matrix[16];

    for (int i=0; i<16 ; i++) matrix[i] = 0;
    matrix[0]  = 1;
    matrix[5]  = 1;
    matrix[10] = 1;
    matrix[15] = 1;
    matrix[12] = x;
    matrix[13] = y;
    matrix[14] = z;
    matXMat(matrix,m,result);
}

}
#endif
