# Project 4: Matrix Multiplication in C
**Name: 周思呈** 
**SID: 12110644** 

## Part 01 - Analysis





## Part 02 - Code
这一部分阐述了一步步优化`matmul_improved` 的过程。
### 访存优化
原本的矩阵乘法外层循环遍历A的行和B的列，最内层遍历A的列和B的行，那么在访问B中元素时就会出现不连续访问的情况，这会导致cache命中率不高，使得数据读写速度变慢。但如果交换j和k的顺序就能巧妙地使B中元素访问连续，提高运算效率。^[# 矩阵乘法&优化方法 - CPU篇 https://zhuanlan.zhihu.com/p/438173915]
![[矩阵乘法.jpeg]]
```c
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
        result->value = (float *) aligned_alloc(128, row * col * sizeof(float));
//#pragma omp parallel for
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
        return result;
    }
}
```

### SIMD指令集
在访存优化的基础上应用SIMD指令集加快计算速度。
SIMD（Single Instruction Multiple Data，单指令多数据流），使用一个控制器控制多个处理单元，同时对一组数据中的每一个数据执行相同的操作。^[# AVX指令集加速矩阵乘法 http://t.csdn.cn/DzZfD]以ARM架构下的NEON指令集为例，主要实现方法为，定义三个能储存四个32位浮点数的寄存器`float32x4_t` ，取出A\[i]\[j]元素并复制四份存在va里，取出B\[j]\[k]以及该地址往后数的四个元素存在vb中，将这二者相乘后加到结果矩阵C\[i]\[k]以及该地址往后数的四个元素vc中。x86架构下的AVX指令集对应寄存器为包含8个float类型的`__m256` ，每次同时对8个元素进行操作。
![[指令集.png]]

```c
Matrix *matmul_improved(const Matrix *matrixLeft, const Matrix *matrixRight)
{
    if (matrixLeft->col != matrixRight->row)
    {
        printf("two matrices can not be multiplied due to illegal size!\n");
        return NULL;
    }
    else
    {
        Matrix *result = (Matrix *) malloc(sizeof(Matrix));
        int row = matrixLeft->row;
        int col = matrixRight->col;
        int mul = matrixLeft->col;
        result->col = col;
        result->row = row;
        result->value = (float *)aligned_alloc(256, row * col * sizeof(float));
        memset(result->value, 0, row * col * sizeof(float));
        float t = 0;
#ifdef WITH_NEON
        float32x4_t va, vb, vc;
        for (int i = 0; i < matrixRight->col; i++) {
            for (int j = 0; j < matrixLeft->row; j++) {
                va = vdupq_n_f32(*(matrixLeft->value + matrixLeft->col * i + j)); //A[i][j]
                for (int k = 0; k < matrixLeft->col; k+=4) {
                    vb = vld1q_f32(matrixRight->value + matrixRight->col * j + k); //B[j][k]
                    vc = vld1q_f32(result->value + result->col * i + k); //C[i][k]
                    vc = vfmaq_f32(vc, vb, va);
                    vst1q_f32(result->value + result->col * i + k, vc);
                }
            }
        }
        return result;
#elifdef WITH_AVX2
        __m256 vecA, vecB, vecC;
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
```




### Strassen算法
![[Pasted image 20221126101012.png]]











## Part 03 - Result & Verifification
### Baseline
进行了指令集的优化之后开始进行第一轮的测试，按照要求测试16x16, 128x128, 1Kx1K, 8Kx8K, 64Kx64K这五种规模的矩阵，于是编写测试代码如下。每种规模的矩阵都测试plain、指令集优化和openblas三种乘法模式的结果，后二者与plain的计算结果进行比较，给出最大误差百分比。
```c
//
// Created by 周思呈 on 2022/11/17.
//
#include "matrix.h"
#include <stdio.h>
#include <sys/time.h>

int main()
{
    struct timeval start,end;
    int size[5] = {16, 128, 1000, 8000, 64000};
    Matrix *matrixx = createMatrix(1000, 1000);
    //热身
    printf("warm up\n");
    Matrix *result = matmul_plain(matrixx, matrixx);
    deleteMatrix(result);
    result = matmul_plain(matrixx, matrixx);
    deleteMatrix(result);

    for (int i = 0; i < 5; ++i) {

        int tempSize = size[i];
        printf("-----------size = %d ---------------\n", tempSize);
        matrixx = createMatrix(tempSize, tempSize);

        printf("this is the plain result\n");
        gettimeofday(&start, NULL);
        result = matmul_plain(matrixx, matrixx);
        gettimeofday(&end, NULL);
        long timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        printf("time=%f\n", timeuse / 1000000.0);
        deleteMatrix(result);

        printf("this is the advanced result\n");
        gettimeofday(&start, NULL);
        Matrix *advanced = matmul_improved(matrixx, matrixx);
        gettimeofday(&end, NULL);
        timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        float mistake = isCorrect(result, advanced);
        printf("time=%f\n", timeuse / 1000000.0);
        printf("the biggest percentage of mistake is %f %%\n", mistake);

        printf("this is the opebBLAS result\n");
        gettimeofday(&start, NULL);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, matrixx->col, matrixx->col, matrixx->col, 1,
                    matrixx->value, matrixx->col, matrixx->value, matrixx->col, 0, advanced->value, matrixx->col);
        gettimeofday(&end, NULL);
        timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
        mistake = isCorrect(result, advanced);
        printf("time=%f\n", timeuse / 1000000.0);
        printf("the biggest percentage of mistake is %f %%\n\n\n", mistake);
        deleteMatrix(advanced);

        deleteMatrix(matrixx);
    }

    return 0;
}
```
但在实际运行过程中发现，8k的朴素矩阵乘法速度非常慢，甚至加了编译优化选项-O3也难以在一两分钟给出计算结果。而16x16和128x128的矩阵在计算速度上差距不是很明显，于是后续测试均选用1k作为矩阵规模。
| 矩阵规模\计算方式  | plain    | openBLAS | neon指令集 |
| ---- | -------- | -------- | ---------- |
| 16   | 0.000023s | 0.000026s | 0.000011s   |
| 128  | 0.011482s | 0.00029s  | 0.005189s   |
| 1000 | 5.902015s | 0.012206s | 2.493775s   |







## Part 04  - Difficulties & Solutions
### 内存管理
在进行优化追求计算速度之前，首先需要保证矩阵乘法计算正确性，因此用规模较小的矩阵作为例子打印其计算结果判断正确性。实现了使用SIMD指令的基本算法之后进行小规模测试，发现前面大部分计算结果都正确，但在最后一行出现了奇怪的错误。
![[内存1.jpg]]
![[内存2.jpg]]
上面一个矩阵是没有使用NEON指令集的普通计算函数的计算结果，下一个是使用之后的。二者对比发现，前面几行计算结果均正确，但最后一行使用指令集之后出现了奇怪的大数。
百思不得其解，让大佬帮忙看代码之后，大佬一针见血地指出：是不是内存没有对齐？看看自己的代码，确实，在`matmul_improved` 里面给结果数组`result->value` 分配内存的时候沿用的是`matmul_plain` 里面的`malloc` ，并不会自动对齐内存。指令集加速的本质是同时操作多个数据，比如simd256就是指同时操作256bit的数据。因此simd256技术要求所操作的数据的首地址是内存对齐32字节。这样的话，原先的malloc就不够用了，因为它不能分配给我们内存对齐于32字节的。^[# 【学习体会】aligned_malloc实现内存对齐 原文链接：https://blog.csdn.net/jin739738709/article/details/122992753] 
用`aligned_alloc` 代替`malloc` 获得连续内存之后获得了正确的计算结果。

### 多核运行计时
在c中计时有多种方法。^[# C语言中常用计时方法总结 http://t.csdn.cn/V9V6z]
一开始本人使用的是最简单的`clock()` ，该函数返回值是硬件滴答数，代码如下。
```cpp
clock_t start,end;
start = clock();
//…calculating…
end = clock();
printf("time=%f\n",(double)end-start)/CLK_TCK);
```
比较`matmul_plain` 和仅添加指令集优化结果时，这个计时方法还能给出相对正确的结果，我们能看到随着矩阵规模的提升，指令集所带来的速度提升越来越大。但是，在使用openMP的情况下，却发现无论是`matmul_plain` 还是`matmul_improved` ，计算耗时都增加了几倍。和同学讨论并查阅资料后发现，`clock()` 计算的是"the CPU time used so far"，所以在打开openMP多个CPU同时运行的时候，他会计算出更长的时间。
换成`gettimeofday()` 之后发现openMP终于给出了令人欣慰的加速结果。
```cpp
struct timeval start,end;
gettimeofday(&start, NULL );
//…calculating…
gettimeofday(&end, NULL );
long timeuse =1000000 * ( end.tv_sec - start.tv_sec ) + end.tv_usec - start.tv_usec;
printf("time=%f\n",timeuse /1000000.0);
```











使用NEON指令集时用到的几个函数
```cpp
float32x4_t res = vdupq_n_f32(0.0f); // 存储的四个float32都初始化为0

float32x4_t res1 = vmlaq_f32(q0, q1, q2); // q0 + q1*q2

float32x4_t q0 = vld1q_f32(d0); // 加载 d0 地址起始的 4 个 float 数据到 q0

```
