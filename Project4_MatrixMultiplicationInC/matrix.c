#include "matrix.h"

void printMatrix(const Matrix *matrix, int precision) //小规模调试
{
    if (!matrix || !matrix->value){
        fprintf(stderr, "matrix to be printed is NULL!\n");
        return;
    }
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
    //但是无法判断地址是否有效
    if (!matrix) {
        fprintf(stderr, "the pointer matrix is NULL!\n");
        return false;
    }
    if (matrix->value) {
        free(matrix->value);
        matrix->value = NULL;
    }
    free(matrix);
    return true;
}

Matrix *createMatrix(const size_t row, const size_t col){
    if (row == 0 || col == 0){
        fprintf(stderr, "rows and/or cols is 0.\n");
        return NULL;
    }
    Matrix *matrix = (Matrix *) malloc(sizeof(Matrix));
    if (matrix == NULL){
        fprintf(stderr, "fail to allocate memory!\n");
        return NULL;
    }
    matrix->col = col;
    matrix->row = row;
//    matrix->value = (float *) aligned_alloc(128, row * col * sizeof(float));
    matrix->value = (float *) malloc(row * col * sizeof(float));
    if (matrix->value == NULL){
        fprintf(stderr, "failed to allocate memory for value!\n");
        free(matrix);
        return NULL;
    }
    memset(matrix->value, 1, col * row * sizeof(float));
    matrix->value[7] = (float)2.2;
    return matrix;
}

float isCorrect(const Matrix *correctMatrix, const Matrix *resultMatrix)
{
    if (!correctMatrix || !resultMatrix){
        fprintf(stderr, "matrix pointer is/are NULL!\n");
        return -1;
    }
    if (!correctMatrix->value || !resultMatrix->value){
        fprintf(stderr, "value pointer is/are NULL!\n");
        return -1;
    }
    if (correctMatrix->col != resultMatrix->col || correctMatrix->row != resultMatrix->col){
        fprintf(stderr, "two matrices have different sizes!\n");
        return -1;
    }
    size_t row = correctMatrix->row;
    size_t col = correctMatrix->col;
    size_t size = row * col;
    float max = 0;
    float original = 0;
    float * p1 = correctMatrix->value;
    float * p2 = resultMatrix->value;
    for (size_t i = 0; i < size; i++){
        float result = *(p1++) - *(p2++);
        if (result > max){
            max = result;
            original = *(p1-1);
        }else if (-result > max){
            max = -result;
            original = *(p1-1);
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

bool matmul_plain(const Matrix *matrixLeft, const Matrix *matrixRight, Matrix * result)
{
    if (!(matrixLeft && matrixRight && result
        && matrixLeft->value && matrixRight->value && result->value)){
        fprintf(stderr, "pointer NULL!\n");
        return 0;
    }
    if (matrixLeft->col != matrixRight->row)
    {
        fprintf(stderr, "two matrices can not be multiplied due to illegal size!\n");
        return 0;
    }
    else
    {
        size_t row = matrixLeft->row;
        size_t col = matrixRight->col;
        size_t mul = matrixLeft->col;
        memset(result->value, 0, row * col * sizeof(float));
#pragma omp parallel for
        for (size_t i = 0; i < row; i++)
        {
            for (size_t k = 0; k < mul; k++)
            {
                for (size_t j = 0; j < col; j++)
                {
                    result->value[i * col + j] += matrixLeft->value[i * mul + k] * matrixRight->value[k * col + j];
                }
            }
        }
        return 1;
    }
}

bool matmul_improved(const Matrix *matrixLeft, const Matrix *matrixRight, Matrix * result)
{
    if (!(matrixLeft && matrixRight && result
          && matrixLeft->value && matrixRight->value && result->value)){
        fprintf(stderr, "pointer NULL!\n");
        return 0;
    }
    if (matrixLeft->col != matrixRight->row)
    {
        fprintf(stderr, "two matrices can not be multiplied due to illegal size!\n");
        return 0;
    }
    else
    {
        size_t M = matrixLeft->row;
        size_t N = matrixRight->col;
        size_t K = matrixLeft->col;
        memset(result->value, 0, M * N * sizeof(float));
        float t = 0;
#ifdef WITH_NEON
        float32x4_t va, vb, vc;
#pragma omp parallel for
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
        return 1;
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
        return 1;
#else
        printf( "NEON and AVX are both not supported\n" );
        return 0;
#endif
    }
}

/*#define UNROLL 4
void dgemm_neon_unroll(size_t n, float *A, float *B, float *C)
{
    for (size_t i = 0; i < n; i += 4*UNROLL) {
        for (size_t j = 0; j < n; j++) {
            float32x4_t c[UNROLL];
            for (int x = 0; x < UNROLL; x++) {
                c[x] = vld1q_f32(C + i + x * 4 + j * n);
            }
            for (size_t k = 0; k < n; k++) {
                float32x4_t b = vdupq_n_f32(B[k + j * n]);
                for (int x = 0; x < UNROLL; x++) {
                    c[x] = vaddq_f32(c[x],
                                     vmulq_f32(vld1q_f32(A + i + k * n + x * 4), b));
                }
            }
            for (int x = 0; x < UNROLL; x++) {
                vst1q_f32(C + i + x * 4 + j * n, c[x]);
            }
        }
    }
}*/

#define UNROLL 4
#define BLOCKSIZE 128

static inline void do_block(size_t n, size_t si, size_t sj, size_t sk,
                            float *A, float *B, float *C)
{
#ifdef WITH_NEON
    for (size_t i = si; i < si + BLOCKSIZE; i += UNROLL*4) {
        for (size_t j = sj; j < sj + BLOCKSIZE; j++) {
            float32x4_t c[UNROLL];
            for (size_t x = 0; x < UNROLL; x++) {
                c[x] = vld1q_f32(C +  i + x * 4 + j * n);
            }
            for (size_t k = sk; k < sk + BLOCKSIZE; k++) {
                float32x4_t b = vdupq_n_f32(B[k + j * n]);
                for (size_t x = 0; x < UNROLL; x++) {
                    c[x] = vaddq_f32(c[x],
                                     vmulq_f32(
                                             vld1q_f32(A + n * k + x * 4 + i), b));
                }
            }

            for (size_t x = 0; x < UNROLL; x++) {
                vst1q_f32(C + i + x * 4 + j * n, c[x]);
            }
        }
    }
#elifdef WITH_AVX2
    for (int i = si; i < si + BLOCKSIZE; i += UNROLL*4) {
		for (int j = sj; j < sj + BLOCKSIZE; j++) {
			__m256d c[UNROLL];
			for (int x = 0; x < UNROLL; x++) {
				c[x] = _mm256_load_pd(C+i+x*4+j*n);
			}
			for (int k = sk; k < sk + BLOCKSIZE; k++) {
				__m256d b = _mm256_broadcast_sd(B+k+j*n);
				for (int x = 0; x < UNROLL; x++) {
					c[x] = _mm256_add_pd(c[x],
						_mm256_mul_pd(
							_mm256_load_pd(A+n*k+x*4+i), b));
				}
			}

			for (int x = 0; x < UNROLL; x++) {
				_mm256_store_pd(C+i+x*4+j*n, c[x]);
			}
		}
	}
#endif
}

void dgemm_neon_unroll_blk(size_t n, float *A, float *B, float *C)
{
#pragma omp parallel for
    for (size_t sj = 0; sj < n; sj += BLOCKSIZE) {
        for (int si = 0; si < n; si += BLOCKSIZE) {
            for (int sk = 0; sk < n; sk += BLOCKSIZE) {
                do_block(n, si, sj, sk, A, B, C);
            }
        }
    }
}

bool neon_unroll(const Matrix *matrixLeft, const Matrix *matrixRight, Matrix *result){
    if (!(matrixLeft && matrixRight && result
          && matrixLeft->value && matrixRight->value && result->value)){
        fprintf(stderr, "pointer NULL!\n");
        return 0;
    }
    if (matrixLeft->col != matrixRight->row)
    {
        fprintf(stderr, "two matrices can not be multiplied due to illegal size!\n");
        return 0;
    }
    size_t M = matrixLeft->row;
    size_t N = matrixRight->col;
    memset(result->value, 0, M * N * sizeof(float));
    dgemm_neon_unroll_blk(M, matrixLeft->value, matrixRight->value, result->value);
    return 1;
}


Matrix *matmul_tile(const Matrix *matrixLeft, const Matrix *matrixRight)
{
    if (matrixLeft->col != matrixRight->row)
    {
        fprintf(stderr, "two matrices can not be multiplied due to illegal size!\n");
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

Matrix* createPlainMatrix(const size_t row, const size_t col){
    if (row == 0 || col == 0){
        fprintf(stderr, "rows and/or cols is 0!\n");
        return NULL;
    }
    Matrix *matrix = (Matrix *) malloc(sizeof(Matrix));
    if (matrix == NULL){
        fprintf(stderr, "failed to allocate memory!\n");
        return NULL;
    }
    matrix->col = col;
    matrix->row = row;
//    matrix->value = (float *) aligned_alloc(128, row * col * sizeof(float));
    matrix->value = (float *) malloc(row * col * sizeof(float));
    if (matrix->value == NULL){
        fprintf(stderr, "failed to allocate memory for value!\n");
        free(matrix);
        return NULL;
    }
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
            vst1q_f32(C->value + indexC + i * C->col + j, vc);
        }
    }
    if (size % 4 != 0){
        for (int j = 0; j < size % 4; j++){
            C->value[indexC + (size-1) * C->col + j] = A->value[indexA + (size-1) * A->col + j] + B->value[indexB + (size-1) * B->col + j];
        }
    }
#elifdef WITH_AVX2
    #pragma omp parallel for
        for (size_t i = 0; i < size; i++)
        {
            for (size_t j = 0; j < size; j += 8)
            {
                _mm256_store_ps((C->value + indexC + i * C->col + j), _mm256_add_ps(_mm256_load_ps(A->value + indexA + i * A->col + j), _mm256_load_ps(B->value + indexB + i * B->col + j)));
            }
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
#elifdef WITH_AVX2
    #pragma omp parallel for
        for (size_t i = 0; i < size; i++)
        {
            for (size_t j = 0; j < size; j += 8)
            {
                _mm256_store_ps((C->value + indexC + i * C->col + j), _mm256_add_ps(_mm256_load_ps(A->value + indexA + i * A->col + j), _mm256_load_ps(B->value + indexB + i * B->col + j)));
            }
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

Matrix *Strassen(size_t indexA, Matrix * A, size_t indexB, Matrix * B, size_t size){
#ifdef WITH_NEON
    if (size <= 256) {
        Matrix *ret = createPlainMatrix(size, size);
        memset(ret->value, 0, ret->col * ret->row * sizeof(float));
        float32x4_t va, vb, vc;
#pragma omp parallel for
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                float t = *(A->value + indexA + A->col * i + j);
                for (int k = 0; k < size; k += 4) {
                    vb = vld1q_f32(B->value + indexB + B->col * j + k); //B[j][k]
                    vc = vld1q_f32(ret->value + ret->col * i + k); //C[i][k]
                    vc = vmlaq_n_f32(vc, vb, t);
                    vst1q_f32(ret->value + ret->col * i + k, vc);
                }
            }
        }
//        dgemm_neon_unroll_blk(size, A->value, B->value, ret->value);
        return ret;
    }
#else
    if (size <= 256){
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
#endif

    //分块
    size_t newSize = (size >> 1);
    size_t a11 = indexA, a12 = indexA + newSize, a21 = a11 + newSize * A->col, a22 = a12 + newSize * A->col;
    size_t b11 = indexB, b12 = indexB + newSize, b21 = b11 + newSize * B->col, b22 = b12 + newSize * B->col;

    //S1 = B12 - B22
    Matrix * S1 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strSubtract(B, b12, B, b22, S1, 0, newSize);
    Matrix * P1 = Strassen(a11, A, 0, S1, newSize);
    deleteMatrix(S1);

    //S3 = A21 + A22
    Matrix * S3 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strAdd(A, a21, A, a22, S3, 0, newSize);
    Matrix * P3 = Strassen(0, S3, b11, B, newSize);
    deleteMatrix(S3);

    //S5 = A11 + A22
    Matrix * S5 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strAdd(A, a11, A, a22, S5, 0, newSize);
    //S6 = B11 + B22
    Matrix * S6 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strAdd(B, b11, B, b22, S6, 0, newSize);
    Matrix * P5 = Strassen(0, S5, 0, S6, newSize);
    deleteMatrix(S5);
    deleteMatrix(S6);

    //S9 = A11 - A21
    Matrix * S9 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strSubtract(A, a11, A, a21, S9, 0, newSize);
    //S10 = B11 + B12
    Matrix * S10 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strAdd(B, b11, B, b12, S10, 0, newSize);
    Matrix * P7 = Strassen(0, S9, 0, S10, newSize);
    deleteMatrix(S9);
    deleteMatrix(S10);

    Matrix * C = createPlainMatrix(size, size);
    size_t c11 = 0, c12 = c11 + newSize, c21 = c11 + newSize * size, c22 = c21 + newSize;
    strAdd(P5, 0, P1, 0, C, c22, newSize);
    strSubtract(C, c22, P3, 0, C, c22, newSize);
    strSubtract(C, c22, P7, 0, C, c22, newSize);
    deleteMatrix(P7);


    //S2 = A11 + A12
    Matrix * S2 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strAdd(A, a11, A, a12, S2, 0, newSize);
    Matrix * P2 = Strassen(0, S2, b22, B, newSize);
    deleteMatrix(S2);
    strAdd(P1, 0, P2, 0, C, c12, newSize);
    deleteMatrix(P1);


    //S4 = B21 - B11
    Matrix * S4 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strSubtract(B, b21, B, b11, S4, 0, newSize);
    Matrix * P4 = Strassen(a22, A, 0, S4, newSize);
    deleteMatrix(S4);
    strAdd(P3, 0, P4, 0, C, c21, newSize);
    deleteMatrix(P3);


    //S7 = A12 - A22
    Matrix * S7 = createPlainMatrix((A->row >> 1), (A->col >> 1));
    strSubtract(A, a12, A, a22, S7, 0, newSize);
    //S8 = B21 + B22
    Matrix * S8 = createPlainMatrix((B->row >> 1), (B->col >> 1));
    strAdd(B, b21, B, b22, S8, 0, newSize);
    Matrix * P6 = Strassen(0, S7, 0, S8, newSize);
    deleteMatrix(S7);
    deleteMatrix(S8);
    strAdd(P5, 0, P4, 0, C, c11, newSize);
    strSubtract(C, c11, P2, 0, C, c11, newSize);
    strAdd(C, c11, P6, 0, C, c11, newSize);

    deleteMatrix(P2);
    deleteMatrix(P4);
    deleteMatrix(P5);
    deleteMatrix(P6);

    return C;

}