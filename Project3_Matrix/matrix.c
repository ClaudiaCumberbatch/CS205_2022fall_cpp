#pragma once

#include "matrix.h"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>

#define delete(p) (deleteMatrix(p), p = NULL)

bool sameSize(const Matrix *matrix1, const Matrix *matrix2)
{
    if ((matrix1->col == matrix2->col && matrix1->row == matrix2->row))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool legalSize(int row, int col)
{
    if (row <= 0 || col <= 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

//判断字符串是否为合法float
bool legal(const char *input)
{
    int *i = (int *)malloc(4);
    *i = 0;
    while (*(input + *i) != 0)
    {
        if ((*(input + *i) < '0' || *(input + *i) > '9') && (*(input + *i) != '.'))
        {
            free(i);
            return false;
        }
        else
        {
            (*i)++;
        }
    }
    free(i);
    return true;
}

Matrix *createMatrix(const int row, const int col, const bool random)
{
    //输入检查
    if (!legalSize(row, col))
    {
        printf("your row and col are too small!\n");
        return NULL;
    }
    else
    {
        Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
        matrix->col = col;
        matrix->row = row;
        char *input = (char *)malloc(sizeof(char) * 64); //用来接收输入的字符串
        // srand((unsigned)time(NULL));
        matrix->value = (float *)malloc(row * col * sizeof(float));
        for (int r = 0; r < row; r++)
        {
            for (int c = 0; c < col; c++)
            {
                if (random)
                {
                    matrix->value[r * col + c] = rand() % 100 + 5;
                }
                else
                {
                    if (c == 0)
                    {
                        printf("please input next %d values of your matrix\n", col);
                    }

                    scanf("%s", input);
                    if (legal(input))
                    {
                        matrix->value[r * col + c] = atof(input);
                    }
                    else
                    {
                        printf("you should input a legal float number\n");
                        c--;
                    }
                }
            }
        }
        free(input);
        return matrix;
    }
}

void printMatrix(const Matrix *matrix, int precision)
{
    if (matrix == NULL)
    {
        printf("this matrix has been deleted!\n");
        return;
    }
    if (precision < 0 || precision > 10)
    {
        printf("the precision is out of range and is automatically set to 2\n");
        precision = 2;
    }

    int col = matrix->col;
    int row = matrix->row;
    printf("the matrix have %d rows and %d cols\n", row, col);
    if (!legalSize(row, col))
    {
        printf("this matrix is empty!\n");
    }
    else
    {
        for (int r = 0; r < row; r++)
        {
            for (int c = 0; c < col; c++)
            {
                printf("%.*lf\t", precision, matrix->value[r * col + c]);
            }
            printf("\n");
        }
    }
}

bool deleteMatrix(Matrix *matrix)
{
    if (matrix != NULL)
    {
        int size = matrix->col * matrix->row;
        for (int i = 0; i < size; i++)
        {
            matrix->value[i] = 0;
        }
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

Matrix *addMatrix(const Matrix *matrix1, const Matrix *matrix2)
{
    if (!(legalSize(matrix1->row, matrix1->col) && legalSize(matrix2->row, matrix2->col)))
    {
        printf("this matrix is empty!\n");
        return NULL;
    }
    if (!sameSize(matrix1, matrix2))
    {
        printf("two matrices have different size which can not be added!\n");
        return NULL;
    }
    // if (!(matrix1->col == matrix2->col && matrix1->row == matrix2->row))
    // {
    //     printf("two matrices have different size which can not be added!\n");
    //     return NULL;
    // }
    // if (matrix1->col <= 0 || matrix1->row <= 0 || matrix2->col < 0 || matrix2->row <= 0)
    // {
    //     printf("this matrix is empty!\n");
    //     return NULL;
    // }
    else
    {
        Matrix *result = (Matrix *)malloc(sizeof(Matrix));
        int row = matrix1->row;
        int col = matrix1->col;
        result->col = col;
        result->row = row;
        int size = col * row;
        result->value = (float *)malloc(row * col * sizeof(float));
        for (int i = 0; i < size; i++)
        {
            if (matrix1->value[i] + matrix2->value[i] > FLT_MAX || matrix1->value[i] + matrix2->value[i] < -FLT_MAX)
            {
                printf("the result is out of range!\n");
                free(result->value);
                free(result);
                return NULL;
            }
            else
            {
                result->value[i] = matrix1->value[i] + matrix2->value[i];
            }
        }
        return result;
    }
}

Matrix *subtractMatrix(const Matrix *matrix1, const Matrix *matrix2)
{
    if (!(legalSize(matrix1->row, matrix1->col) && legalSize(matrix2->row, matrix2->col)))
    {
        printf("this matrix is empty!\n");
        return NULL;
    }
    if (!sameSize(matrix1, matrix2))
    {
        printf("two matrices have different size which can not be subtracted!\n");
        return NULL;
    }
    else
    {
        Matrix *result = (Matrix *)malloc(sizeof(Matrix));
        int row = matrix1->row;
        int col = matrix1->col;
        result->col = col;
        result->row = row;
        int size = col * row;
        result->value = (float *)malloc(row * col * sizeof(float));
        for (int i = 0; i < size; i++)
        {
            if (matrix1->value[i] - matrix2->value[i] > FLT_MAX || matrix1->value[i] - matrix2->value[i] < -FLT_MAX)
            {
                printf("the result is out of range!\n");
                free(result->value);
                free(result);
                return NULL;
            }
            else
            {
                result->value[i] = matrix1->value[i] - matrix2->value[i];
            }
        }
        return result;
    }
}

Matrix *operateScalar(const Matrix *matrix, const float scalar, const int opr)
{
    if (!(legalSize(matrix->row, matrix->col)))
    {
        printf("this matrix is empty!\n");
        return NULL;
    }
    // if (matrix->col <= 0 || matrix->row <= 0)
    // {
    //     printf("your row and col are too small!\n");
    //     return NULL;
    // }

    Matrix *result = (Matrix *)malloc(sizeof(Matrix));
    int row = matrix->row;
    int col = matrix->col;
    result->col = col;
    result->row = row;
    int size = col * row;
    result->value = (float *)malloc(row * col * sizeof(float));
    for (int i = 0; i < size; i++)
    {
        if (((opr == ADD) && (matrix->value[i] + scalar > FLT_MAX || matrix->value[i] + scalar < -FLT_MAX)) || ((opr == MINUS) && (matrix->value[i] - scalar > FLT_MAX || matrix->value[i] - scalar < -FLT_MAX)) || ((opr == MULTIPLY) && (matrix->value[i] * scalar > FLT_MAX || matrix->value[i] * scalar < -FLT_MAX)))
        {
            printf("the result is out of range!\n");
            free(result->value);
            free(result);
            return NULL;
        }
        else
        {
            if (opr == ADD)
            {
                result->value[i] = matrix->value[i] + scalar;
            }
            else if (opr == MINUS)
            {
                result->value[i] = matrix->value[i] - scalar;
            }
            else if (opr == MULTIPLY)
            {
                result->value[i] = matrix->value[i] * scalar;
            }
            else
            {
                printf("wrong operator!\n");
                free(result->value);
                free(result);
                return NULL;
            }
        }
    }
    return result;
}

Matrix *multiplyMatrix(const Matrix *matrixLeft, const Matrix *matrixRight)
{
    if (matrixLeft->col != matrixRight->row)
    {
        printf("two matrices can not be multiplied due to illegal size!\n");
        return NULL;
    }
    else
    {
        Matrix *result = (Matrix *)malloc(sizeof(Matrix));
        int row = matrixLeft->row;
        int col = matrixRight->col;
        int mul = matrixLeft->col;
        result->col = col;
        result->row = row;
        result->value = (float *)malloc(row * col * sizeof(float));
        float t = 0;
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                t = 0;
                for (int k = 0; k < mul; k++)
                {
                    if (t + matrixLeft->value[i * mul + k] * matrixRight->value[k * col + j] > FLT_MAX || t + matrixLeft->value[i * mul + k] * matrixRight->value[k * col + j] < -FLT_MAX)
                    {
                        printf("the result is out of range!\n");
                        free(result->value);
                        free(result);
                        return NULL;
                    }
                    else
                    {
                        t += matrixLeft->value[i * mul + k] * matrixRight->value[k * col + j];
                    }
                }
                result->value[i * col + j] = t;
            }
        }
        return result;
    }
}

float extremeValue(const Matrix *matrix, bool max)
{
    if (matrix == NULL)
    {
        printf("this matrix has been deleted!\n");
        return 0;
    }
    int col = matrix->col;
    int row = matrix->row;
    if (!legalSize(row, col))
    {
        printf("this matrix is empty!\n");
        return 0;
    }
    else
    {
        int size = col * row;
        float result = matrix->value[0];
        for (size_t i = 1; i < size; i++)
        {
            if (max)
            {
                if (matrix->value[i] > result)
                {
                    result = matrix->value[i];
                }
            }
            else
            {
                if (matrix->value[i] < result)
                {
                    result = matrix->value[i];
                }
            }
        }

        // for (int r = 0; r < row; r++)
        // {
        //     for (int c = 0; c < col; c++)
        //     {
        //         if (max)
        //         {
        //             if (matrix->value[r * col + c] > result)
        //             {
        //                 result = matrix->value[r * col + c];
        //             }
        //         }
        //         else
        //         {
        //             if (matrix->value[r * col + c] < result)
        //             {
        //                 result = matrix->value[r * col + c];
        //             }
        //         }
        //     }
        // }
        return result;
    }
}

Matrix *transpose(const Matrix *matrix)
{
    if (matrix == NULL)
    {
        printf("this matrix has been deleted!\n");
        return NULL;
    }
    int col = matrix->col;
    int row = matrix->row;
    if (row <= 0 || col <= 0)
    {
        printf("this matrix is empty!\n");
        return NULL;
    }
    else
    {
        Matrix *result = (Matrix *)malloc(sizeof(Matrix));
        result->col = row;
        result->row = col;
        result->value = (float *)malloc(row * col * sizeof(float));
        for (int i = 0; i < result->row; i++)
        {
            for (int j = 0; j < result->col; j++)
            {
                result->value[i * result->col + j] = matrix->value[j * col + i];
            }
        }
        return result;
    }
}

Matrix *identityMatrix(const int sideLength)
{
    //输入检查
    if (sideLength <= 0)
    {
        printf("your row and col are too small!\n");
        return NULL;
    }
    else
    {
        Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
        matrix->col = sideLength;
        matrix->row = sideLength;
        matrix->value = (float *)malloc(sideLength * sideLength * sizeof(float));
        for (int r = 0; r < sideLength; r++)
        {
            for (int c = 0; c < sideLength; c++)
            {
                if (r == c)
                {
                    matrix->value[r * sideLength + c] = 1;
                }
                else
                {
                    matrix->value[r * sideLength + c] = 0;
                }
            }
        }
        return matrix;
    }
}

bool set(Matrix *matrix, const int r, const int c, const float newValue)
{
    if (matrix == NULL)
    {
        printf("this matrix has been deleted!\n");
        return false;
    }
    if (matrix->col < c - 1 || matrix->row < r - 1)
    {
        printf("this place is illegal!\n");
        return false;
    }
    matrix->value[(r - 1) * matrix->col + (c - 1)] = newValue;
    return true;
}

Matrix *copy(const Matrix *matrix)
{
    if (matrix == NULL)
    {
        printf("this matrix has been deleted!\n");
        return NULL;
    }
    else
    {
        Matrix *result = (Matrix *)malloc(sizeof(Matrix));
        int row = matrix->row;
        int col = matrix->col;
        result->col = col;
        result->row = row;
        int size = col * row;
        result->value = (float *)malloc(row * col * sizeof(float));
        for (int i = 0; i < size; i++)
        {
            result->value[i] = matrix->value[i];
        }
        return result;
    }
}

Matrix *ones(const int row, const int col)
{
    if (row <= 0 || col <= 0)
    {
        printf("your row and col are too small!\n");
        return NULL;
    }
    Matrix *result = (Matrix *)malloc(sizeof(Matrix));
    result->col = col;
    result->row = row;
    int size = col * row;
    result->value = (float *)malloc(row * col * sizeof(float));
    for (int i = 0; i < size; i++)
    {
        result->value[i] = 1;
    }
    return result;
}

Matrix *zeros(const int row, const int col)
{
    if (row <= 0 || col <= 0)
    {
        printf("your row and col are too small!\n");
        return NULL;
    }
    Matrix *result = (Matrix *)malloc(sizeof(Matrix));
    result->col = col;
    result->row = row;
    int size = col * row;
    result->value = (float *)malloc(row * col * sizeof(float));
    for (int i = 0; i < size; i++)
    {
        result->value[i] = 0;
    }
    return result;
}