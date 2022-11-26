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
    printf("warm up\n");
    Matrix *result = matmul_plain(matrixx, matrixx);
//    deleteMatrix(result);
    result = matmul_plain(matrixx, matrixx);
//    deleteMatrix(result);

//    for (int i = 0; i < 5; ++i) {

        int tempSize = 1024;
        printf("-----------size = %d ---------------\n", tempSize);
        matrixx = createMatrix(tempSize, tempSize);
//        printMatrix(matrixx,0);

        printf("this is the plain result\n");
        gettimeofday(&start, NULL);
        result = matmul_plain(matrixx, matrixx);
        gettimeofday(&end, NULL);
        long timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        printf("time=%f\n", timeuse / 1000000.0);
//        deleteMatrix(result);
//        printMatrix(result, 0);

        printf("this is the advanced result\n");
        gettimeofday(&start, NULL);
//        Matrix *advanced = matmul_improved(matrixx, matrixx);
        Matrix *advanced = Strassen(0, matrixx, 0, matrixx, matrixx->col);
        gettimeofday(&end, NULL);
        timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        float mistake = isCorrect(result, advanced);
        printf("time=%f\n", timeuse / 1000000.0);
        printf("the biggest percentage of mistake is %f %%\n", mistake);
//        printMatrix(advanced, 0);

        printf("this is the opebBLAS result\n");
        gettimeofday(&start, NULL);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, matrixx->col, matrixx->col, matrixx->col, 1,
                    matrixx->value, matrixx->col, matrixx->value, matrixx->col, 0, advanced->value, matrixx->col);
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