#pragma once

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>

#define ADD 1
#define MINUS 2
#define MULTIPLY 3

typedef struct _Matrix
{
    int row;
    int col;
    float *value;
} Matrix;

bool legal(const char *input);

Matrix *createMatrix(const int row, const int col, const bool random);

void printMatrix(const Matrix *matrix, int precision);

bool deleteMatrix(Matrix *matrix);

Matrix *addMatrix(const Matrix *matrix1, const Matrix *matrix2);

Matrix *subtractMatrix(const Matrix *matrix1, const Matrix *matrix2);

Matrix *operateScalar(const Matrix *matrix, const float scalar, const int opr);

Matrix *multiplyMatrix(const Matrix *matrixLeft, const Matrix *matrixRight);

float extremeValue(const Matrix *matrix, bool max);

Matrix *transpose(const Matrix *matrix);

Matrix *identityMatrix(const int sideLength);

bool set(Matrix *matrix, const int r, const int c, const float newValue);

Matrix *copy(const Matrix *matrix);

Matrix *ones(const int row, const int col);

Matrix *zeros(const int row, const int col);