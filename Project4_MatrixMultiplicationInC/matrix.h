#pragma once

#define WITH_NEON
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include </opt/homebrew/Cellar/openblas/0.3.21/include/cblas.h>
#include <string.h>
#include <math.h>

#define _OMP_THREAD_ 4

#ifdef WITH_NEON
#include <arm_neon.h>
#endif

#ifdef WITH_AVX2
#include <immintrin.h>
#endif

#ifdef _OPENMP
#include </opt/homebrew/opt/libomp/include/omp.h>
#endif

typedef struct _Matrix
{
    size_t row;
    size_t col;
    float *value;
} Matrix;

void printMatrix(const Matrix *matrix, int precision);
Matrix *createMatrix(const size_t row, const size_t col);
bool deleteMatrix(Matrix *matrix);
float isCorrect(const Matrix *correctMatrix, const Matrix *resultMatrix);
Matrix *matmul_plain(const Matrix *matrixLeft, const Matrix *matrixRight);
Matrix *matmul_improved(const Matrix *matrixLeft, const Matrix *matrixRight);
Matrix *matmul_tile(const Matrix *matrixLeft, const Matrix *matrixRight);

Matrix* createPlainMatrix(size_t row, size_t col);
void strAdd(Matrix * A, size_t indexA, Matrix * B, size_t indexB, Matrix * C, size_t indexC, size_t size);
void strSubtract(Matrix * A, size_t indexA, Matrix * B, size_t indexB, Matrix * C, size_t indexC, size_t size);
Matrix *Strassen(size_t indexA, Matrix * A, size_t indexB, Matrix * B, size_t size);

