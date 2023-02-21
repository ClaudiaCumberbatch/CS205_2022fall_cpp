#include "matrix.c"
// #include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <time.h>
clock_t start, end;

int main()
{
    start = clock();

    // Matrix *matrix = createMatrix(3, 4, 1);
    // printMatrix(matrix, 0);
    // delete (matrix);
    // printMatrix(matrix, 0);

    // Matrix *matrix = createMatrix(3, 4, 1);
    // printMatrix(matrix, 0);
    // Matrix *matrix2 = copy(matrix);
    // printMatrix(matrix2, 0);

    // bool flag = deletee(matrix);
    // matrix = NULL;
    // printf("1st time %s\n", flag?"success":"no");
    // printMatrix(matrix);

    // flag = deletee(matrix);
    // printf("2nd time %s\n", flag?"success":"no");

    Matrix *matrix1 = createMatrix(10, 10, 1);
    Matrix *matrix2 = createMatrix(10, 10, 1);
    Matrix *matrix3 = createMatrix(1, 2, 1);
    printMatrix(matrix1, 0);
    printf("\n");
    printMatrix(matrix2, 0);
    printf("\n");
    printMatrix(matrix3, 0);
    printf("\n");
    Matrix *addresult1 = addMatrix(matrix1, matrix2);
    printMatrix(addresult1, 0);
    printf("\n");
    Matrix *addresult2 = addMatrix(matrix1, matrix3);
    printMatrix(addresult2, 0);
    printf("\n");
    Matrix *minusresult = subtractMatrix(matrix1, matrix2);
    printMatrix(minusresult, 0);

    // float scalar = 2.333;
    // Matrix *addresult = operateScalar(matrix, scalar, ADD);
    // printMatrix(addresult,1);
    // Matrix *minusresult = operateScalar(matrix, scalar, MINUS);
    // printMatrix(minusresult,1);
    // Matrix *mulresult = operateScalar(matrix, scalar, MULTIPLY);
    // printMatrix(mulresult,1);

    // Matrix *matrix1 = createMatrix(4, 3, 1);
    // Matrix *matrix2 = createMatrix(3, 5, 1);
    // printMatrix(matrix1, 0);
    // printf("\n");
    // printMatrix(matrix2,0);
    // printf("\n");
    // printMatrix(multiplyMatrix(matrix1, matrix2),0);
    // printf("\n");
    // printMatrix(multiplyMatrix(matrix2, matrix1),0);
    // printf("\n");

    // printf("%f\n", extremeValue(matrix1, 1));
    // printf("%f\n", extremeValue(matrix1, 0));
    // printMatrix(transpose(matrix1), 0);

    // printMatrix(identityMatrix(6), 0);
    // printMatrix(ones(3, 4), 0);
    // printMatrix(zeros(4, 5), 0);

    // Matrix *matrix2 = createMatrix(4, 3, 1);
    // printMatrix(matrix2, 1);
    // set(matrix2, 1, 2, 6.6);
    // printMatrix(matrix2, 1);

    // printMatrix(ones(3,4),2);

    end = clock(); //程序结束用时
    // double endtime = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Total time: %f ms\n", endtime * 1000); // ms为单位

    return 0;
}