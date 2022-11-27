//
// Created by 周思呈 on 2022/11/17.
//
#include "matrix.h"
#include <stdio.h>
#include <sys/time.h>

int main()
{
    struct timeval start,end;
//    int size[5] = {16, 128, 1000, 8000, 64000};
    Matrix *matrixx = createMatrix(1000, 1000);
    //热身
    Matrix *result = createPlainMatrix(1000, 1000);
    printf("warm up\n");
    if (matmul_plain(matrixx, matrixx, result)){
        printf("warm up 2 done\n");
    }
    deleteMatrix(result);
    result = NULL;
    result = createPlainMatrix(1000, 1000);
    if (matmul_plain(matrixx, matrixx, result)){
        printf("warm up 1 done\n");
    }
    deleteMatrix(result);
    result = NULL;
    deleteMatrix(matrixx);
    matrixx = NULL;

//    for (int i = 0; i < 5; ++i) {

        int tempSize = 1024;
        printf("-----------size = %d ---------------\n", tempSize);
        matrixx = createMatrix(tempSize, tempSize);
//        printMatrix(matrixx,0);

        printf("this is the plain result\n");
        result = createPlainMatrix(tempSize, tempSize);
        gettimeofday(&start, NULL);
        matmul_plain(matrixx, matrixx, result);
        gettimeofday(&end, NULL);
        long timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        printf("time=%f\n", timeuse / 1000000.0);
//        printMatrix(result, 0);

        printf("this is the advanced result\n");
        Matrix *advanced = createPlainMatrix(tempSize, tempSize);
        gettimeofday(&start, NULL);
        neon_unroll(matrixx, matrixx, advanced);
        gettimeofday(&end, NULL);
        timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        float mistake = isCorrect(result, advanced);
        printf("time=%f\n", timeuse / 1000000.0);
        printf("the biggest percentage of mistake is %f %%\n", mistake);
        deleteMatrix(advanced);
//        printMatrix(advanced, 0);

        printf("this is the opebBLAS result\n");
        advanced = createPlainMatrix(tempSize, tempSize);
        gettimeofday(&start, NULL);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, matrixx->col, matrixx->col, matrixx->col, 1,
                    matrixx->value, matrixx->col, matrixx->value, matrixx->col, 0, advanced->value, matrixx->col);
        gettimeofday(&end, NULL);
        timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        mistake = isCorrect(result, advanced);
        printf("time=%f\n", timeuse / 1000000.0);
        printf("the biggest percentage of mistake is %f %%\n", mistake);
        deleteMatrix(advanced);

        printf("this is the strassen result\n");
        gettimeofday(&start, NULL);
        advanced = Strassen(0, matrixx, 0, matrixx, matrixx->col);
        gettimeofday(&end, NULL);
        timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        mistake = isCorrect(result, advanced);
        printf("time=%f\n", timeuse / 1000000.0);
        printf("the biggest percentage of mistake is %f %%\n", mistake);

        deleteMatrix(advanced);
        deleteMatrix(matrixx);
//    }

    return 0;
}