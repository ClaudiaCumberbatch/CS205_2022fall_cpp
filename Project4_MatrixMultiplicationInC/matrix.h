#ifndef _MATRIX_H
#define _MATRIX_H

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
bool matmul_plain(const Matrix *matrixLeft, const Matrix *matrixRight, Matrix * result);
bool matmul_improved(const Matrix *matrixLeft, const Matrix *matrixRight, Matrix * result);
Matrix *matmul_tile(const Matrix *matrixLeft, const Matrix *matrixRight);

void dgemm_neon_unroll(size_t n, float *A, float *B, float *C);

static inline void do_block(size_t n, size_t si, size_t sj, size_t sk, float *A, float *B, float *C);
void dgemm_neon_unroll_blk(size_t n, float *A, float *B, float *C);
bool neon_unroll(const Matrix *matrixLeft, const Matrix *matrixRight, Matrix *result);

Matrix* createPlainMatrix(size_t row, size_t col);
void strAdd(Matrix * A, size_t indexA, Matrix * B, size_t indexB, Matrix * C, size_t indexC, size_t size);
void strSubtract(Matrix * A, size_t indexA, Matrix * B, size_t indexB, Matrix * C, size_t indexC, size_t size);
Matrix *Strassen(size_t indexA, Matrix * A, size_t indexB, Matrix * B, size_t size);

#endif