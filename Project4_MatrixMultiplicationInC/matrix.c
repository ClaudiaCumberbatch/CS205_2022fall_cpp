#include "matrix.h"

void printMatrix(const Matrix *matrix, int precision)
{
    size_t col = matrix->col;
    size_t row = matrix->row;
    printf("the matrix have %d rows and %d cols\n", row, col);

    for (size_t r = 0; r < row; r++)
    {
        for (size_t c = 0; c < col; c++)
        {
            printf("%.*lf\t", precision, matrix->value[r * col + c]);
        }
        printf("\n");
    }
}

bool deleteMatrix(Matrix *matrix)
{
    if (matrix != NULL)
    {
        size_t size = matrix->col * matrix->row;
        memset(matrix->value, 0, size * sizeof(float));
        free(matrix->value);
        matrix->value = NULL;
        matrix->col = 0;
        matrix->row = 0;
        free(matrix);
        // printf("please remember to set this pointer to NULL\n");
        return true;
    }
    else
    {
        return false;
    }
}

Matrix *createMatrix(const size_t row, const size_t col){
    Matrix *matrix = (Matrix *) malloc(sizeof(Matrix));
    matrix->col = col;
    matrix->row = row;
//    matrix->value = (float *) aligned_alloc(128, row * col * sizeof(float));
    matrix->value = (float *) malloc(row * col * sizeof(float));
//    memset(matrix->value, (float)3, matrix->col * matrix->row);
//    matrix->value[6] = 2;
    for (int r = 0; r < row; r++)
    {
        for (int c = 0; c < col; c++)
        {
            if (c == 1){
                matrix->value[r * col + c] = 2;
            }
            else{
                matrix->value[r * col + c] = 1;
            }
//            matrix->value[r * col + c] = r * col + c;
        }
    }
    return matrix;
}

float isCorrect(const Matrix *correctMatrix, const Matrix *resultMatrix)
{
    size_t row = correctMatrix->row;
    size_t col = correctMatrix->col;
    float max = 0;
    float original = 0;
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            if (*(correctMatrix->value + i * col + j) - *(resultMatrix->value + i * col + j) > max)
            {
                max = *(correctMatrix->value + i * col + j) - *(resultMatrix->value + i * col + j);
                original = *(correctMatrix->value + i * col + j);
            }else if(-*(correctMatrix->value + i * col + j) + *(resultMatrix->value + i * col + j) > max)
            {
                max = -*(correctMatrix->value + i * col + j) + *(resultMatrix->value + i * col + j);
                original = *(resultMatrix->value + i * col + j);
            }
        }
    }
    if (original == 0)
    {
        return 0;
    }
    else
    {
        printf("the correct result is %f while the biggest difference is %f\n", original, max);
        return max*(float)100/original;
    }
}

Matrix *matmul_plain(const Matrix *matrixLeft, const Matrix *matrixRight)
{
    if (matrixLeft->col != matrixRight->row)
    {
        printf("two matrices can not be multiplied due to illegal size!\n");
        return NULL;
    }
    else
    {
        Matrix *result = (Matrix *) malloc(sizeof(Matrix));
        size_t row = matrixLeft->row;
        size_t col = matrixRight->col;
        size_t mul = matrixLeft->col;
        result->col = col;
        result->row = row;
//        result->value = (float *) aligned_alloc(128, row * col * sizeof(float));
        result->value = (float *) malloc( row * col * sizeof(float));
#pragma omp parallel for
        for (size_t i = 0; i < row; i++)
        {
            for (size_t k = 0; k < mul; k++)
            {
                for (size_t j = 0; j < col; j++)
                {
                    /*if (matrixLeft->col == 10){
                        printf("i = %d, j = %d, k = %d\n", i, j, k);
                    }*/
                    result->value[i * col + j] += matrixLeft->value[i * mul + k] * matrixRight->value[k * col + j];
                }
            }
        }
        return result;
    }
}

Matrix *matmul_improved(const Matrix *matrixLeft, const Matrix *matrixRight)
{
    if (matrixLeft->col != matrixRight->row)
    {
        printf("two matrices can not be multiplied due to illegal size!\n");
        printf("left col is %d while right row is %d\n", matrixLeft->col, matrixRight->col);
        return NULL;
    }
    else
    {
        Matrix *result = (Matrix *) malloc(sizeof(Matrix));
        size_t M = matrixLeft->row;
        size_t N = matrixRight->col;
        size_t K = matrixLeft->col;
        result->col = N;
        result->row = M;
        result->value = (float *)aligned_alloc(256, M * N * sizeof(float));
        memset(result->value, 0, M * N * sizeof(float));
        float t = 0;
#ifdef WITH_NEON
        float32x4_t va, vb, vc;
//#pragma omp parallel for
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                float t = *(matrixLeft->value + M * i + j);
                for (int k = 0; k < K; k+=4) {
                    vb = vld1q_f32(matrixRight->value + N * j + k); //B[j][k]
                    vc = vld1q_f32(result->value + N * i + k); //C[i][k]
                    vc = vmlaq_n_f32(vc, vb, t);
                    vst1q_f32(result->value + N * i + k, vc);
                }
            }
        }

//        for (int i = 0; i < matrixRight->col; i++) {
//            for (int j = 0; j < matrixLeft->row; j++) {
//                va = vdupq_n_f32(*(matrixLeft->value + matrixLeft->col * i + j)); //A[i][j]
//                for (int k = 0; k < matrixLeft->col; k+=4) {
//                    vb = vld1q_f32(matrixRight->value + matrixRight->col * j + k); //B[j][k]
//                    vc = vld1q_f32(result->value + result->col * i + k); //C[i][k]
//                    vc = vfmaq_f32(vc, vb, va);
//                    vst1q_f32(result->value + result->col * i + k, vc);
//                }
//            }
//        }
        return result;
#elifdef WITH_AVX2
        __m256 vecA, vecB, vecC;
#pragma omp parallel for
        for (int i = 0; i < matrixRight->col; i++) {
            for (int j = 0; j < matrixLeft->row; j++) {
                vecA = _mm256_set1_ps(*(matrixLeft->value + matrixLeft->col * i + j)); //A[i][j]
                for (int k = 0; k < matrixLeft->col; k+=8) {
                    vecB = _mm256_loadu_ps(matrixRight->value + matrixRight->col * j + k); //B[j][k]
                    vecC = _mm256_loadu_ps(result->value + result->col * i + k); //C[i][k]
                    vecC = _mm256_fmadd_ps(va, vb, vc);
                    _mm256_storeu_ps(result->value + result->col * i + k, vecC);
                }
            }
        }
        return result;
#else
        printf( "NEON and AVX are both not supported\n" );
        return NULL;
#endif
    }
}

Matrix *matmul_tile(const Matrix *matrixLeft, const Matrix *matrixRight)
{
    if (matrixLeft->col != matrixRight->row)
    {
        printf("two matrices can not be multiplied due to illegal size!\n");
        return NULL;
    }
    else
    {
        Matrix *result = (Matrix *) malloc(sizeof(Matrix));
        size_t M = matrixLeft->row;
        size_t N = matrixRight->col;
        size_t K = matrixLeft->col;
        result->col = M;
        result->row = N;
        result->value = (float *) aligned_alloc(128, M * N * sizeof(float));
        size_t iTile = 125, jTile = 125, kTile = 125;
        size_t iOuterBound = M/iTile, jOuterBound = N/jTile, kOuterBound = K/kTile;
//#pragma omp parallel for
        for (size_t i_outer = 0; i_outer < iOuterBound; i_outer++) {
            for (int j_outer = 0; j_outer < jOuterBound; j_outer++) {
                for (int k_outer = 0; k_outer < kOuterBound; k_outer++) {
                    for (int i_inner = 0; i_inner < iTile; i_inner++) {
                        for (int k_inner = 0; k_inner < kTile; k_inner++) {
                            for (int j_inner = 0; j_inner < jTile; j_inner++) {
                                result->value[(i_outer * iTile + i_inner) * N +
                                  (j_outer * jTile + j_inner)] +=
                                        matrixLeft->value[(i_outer * iTile + i_inner) * K +
                                          (k_outer * kTile + k_inner)] *
                                        matrixRight->value[(k_outer * kTile + k_inner) * N +
                                          (j_outer * jTile + j_inner)];
                            }
                        }
                    }
                }
            }
        }
        return result;
    }
}

Matrix* createPlainMatrix(size_t row, size_t col){
    Matrix *matrix = (Matrix *) malloc(sizeof(Matrix));
    matrix->col = col;
    matrix->row = row;
//    matrix->value = (float *) aligned_alloc(128, row * col * sizeof(float));
    matrix->value = (float *) malloc(row * col * sizeof(float));
//    printf("%x\n", matrix->value);
//    printf("%f\n", matrix->value[0]);
    memset(matrix->value, 0, col * row * sizeof(float));

    return matrix;
}

void strAdd(Matrix * A, size_t indexA, Matrix * B, size_t indexB, Matrix * C, size_t indexC, size_t size){
#ifdef WITH_NEON
    float32x4_t va, vb, vc;
#pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; j+=4){
            va = vld1q_f32(A->value + indexA + i * A->col + j);
            vb = vld1q_f32(B->value + indexB + i * B->col + j);
            vc = vaddq_f32(va, vb);
//#pragma omp critical
            vst1q_f32(C->value + indexC + i * C->col + j, vc);
        }
    }
    if (size % 4 != 0){
        for (int j = 0; j < size % 4; j++){
            C->value[indexC + (size-1) * C->col + j] = A->value[indexA + (size-1) * A->col + j] + B->value[indexB + (size-1) * B->col + j];
        }
    }
#else
    printf( "NEON and AVX are both not supported\n" );
        return NULL;
#endif
}

void strSubtract(Matrix * A, size_t indexA, Matrix * B, size_t indexB, Matrix * C, size_t indexC, size_t size){
#ifdef WITH_NEON
    float32x4_t va, vb, vc;
#pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; j+=4){
            va = vld1q_f32(A->value + indexA + i * A->col + j);
            vb = vld1q_f32(B->value + indexB + i * B->col + j);
            vc = vsubq_f32(va, vb);
//#pragma omp critical
            vst1q_f32(C->value + indexC + i * C->col + j, vc);
        }
    }
    if (size % 4 != 0){
        for (int j = 0; j < size % 4; j++){
            C->value[indexC + (size-1) * C->col + j] = A->value[indexA + (size-1) * A->col + j] - B->value[indexB + (size-1) * B->col + j];
        }
    }
#else
    printf( "NEON and AVX are both not supported\n" );
        return NULL;
#endif
}

Matrix *Strassen(size_t indexA, Matrix * A, size_t indexB, Matrix * B, size_t size){
    if (size <= 128){
        Matrix * ret = createPlainMatrix(size, size);
        memset(ret->value, 0, ret->col * ret->row * sizeof(float));
#pragma omp parallel for
        for (size_t i = 0; i < size; i++)
        {
            for (size_t k = 0; k < size; k++)
            {
                for (size_t j = 0; j < size; j++)
                {
                    ret->value[i * ret->col + k] += A->value[indexA + i * A->col + j] * B->value[indexB + j * B->col + k];
                }
            }
        }
        return ret;
    }

    //分块
    size_t newSize = (size >> 1);
    size_t a11 = indexA, a12 = indexA + newSize, a21 = a11 + newSize * A->col, a22 = a12 + newSize * A->col;
    size_t b11 = indexB, b12 = indexB + newSize, b21 = b11 + newSize * B->col, b22 = b12 + newSize * B->col;

    //S1 = B12 - B22
    Matrix * S1 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strSubtract(B, b12, B, b22, S1, 0, newSize);

    //S2 = A11 + A12
    Matrix * S2 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strAdd(A, a11, A, a12, S2, 0, newSize);

    //S3 = A21 + A22
    Matrix * S3 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strAdd(A, a21, A, a22, S3, 0, newSize);

    //S4 = B21 - B11
    Matrix * S4 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strSubtract(B, b21, B, b11, S4, 0, newSize);

    //S5 = A11 + A22
    Matrix * S5 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strAdd(A, a11, A, a22, S5, 0, newSize);

    //S6 = B11 + B22
    Matrix * S6 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strAdd(B, b11, B, b22, S6, 0, newSize);

    //S7 = A12 - A22
    Matrix * S7 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strSubtract(A, a12, A, a22, S7, 0, newSize);

    //S8 = B21 + B22
    Matrix * S8 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strAdd(B, b21, B, b22, S8, 0, newSize);

    //S9 = A11 - A21
    Matrix * S9 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strSubtract(A, a11, A, a21, S9, 0, newSize);

    //S10 = B11 + B12
    Matrix * S10 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strAdd(B, b11, B, b12, S10, 0, newSize);

    //P1 = A11 * S1
    Matrix * P1 = Strassen(a11, A, 0, S1, newSize);
    Matrix * P2 = Strassen(0, S2, b22, B, newSize);
    Matrix * P3 = Strassen(0, S3, b11, B, newSize);
    Matrix * P4 = Strassen(a22, A, 0, S4, newSize);
    Matrix * P5 = Strassen(0, S5, 0, S6, newSize);
    Matrix * P6 = Strassen(0, S7, 0, S8, newSize);
    Matrix * P7 = Strassen(0, S9, 0, S10, newSize);

    Matrix * C = createPlainMatrix(size, size);
    size_t c11 = 0, c12 = c11 + newSize, c21 = c11 + newSize * size, c22 = c21 + newSize;

    strAdd(P5, 0, P4, 0, C, c11, newSize);
    strSubtract(C, c11, P2, 0, C, c11, newSize);
    strAdd(C, c11, P6, 0, C, c11, newSize);

    strAdd(P1, 0, P2, 0, C, c12, newSize);

    strAdd(P3, 0, P4, 0, C, c21, newSize);

    strAdd(P5, 0, P1, 0, C, c22, newSize);
    strSubtract(C, c22, P3, 0, C, c22, newSize);
    strSubtract(C, c22, P7, 0, C, c22, newSize);

    deleteMatrix(S1);
    deleteMatrix(S2);
    deleteMatrix(S3);
    deleteMatrix(S4);
    deleteMatrix(S5);
    deleteMatrix(S6);
    deleteMatrix(S7);
    deleteMatrix(S8);
    deleteMatrix(S9);
    deleteMatrix(S10);
    deleteMatrix(P1);
    deleteMatrix(P2);
    deleteMatrix(P3);
    deleteMatrix(P4);
    deleteMatrix(P5);
    deleteMatrix(P6);
    deleteMatrix(P7);

    return C;

}