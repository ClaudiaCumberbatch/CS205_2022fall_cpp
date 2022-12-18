//
// Created by 周思呈 on 2022/12/12.
//
#include "Matrix.h"
using namespace std;

int main()
{
    /*int* ip = new int[12 * sizeof(int)]{1,1,1,1,
                                        2,2,2,2,
                                        3,3,3,3};
    double* dp = new double[12 * sizeof(double)]{1.1,2.2,3.3,4.4,
                                                 5.5,6.6,7.7,8.8,
                                                 9.9,1.1,2.2,3.3};
    Matrix matrix(Matrix::TYPE4I, 3, 4, ip,1);
    int i = 0;
    matrix.getElement(6,3,1,i);
    matrix.setElement(3,6,1,i);


    cout << matrix+1 << endl;
    cout << matrix-1 << endl;
    Matrix matrix2(Matrix::TYPE64F, 4, 3, dp,1);
    Matrix mt = matrix * matrix2;
    cout << matrix << endl;
    cout << matrix2 << endl;
    cout << mt << endl;
    cout << mt.TYPE << endl;
    visit([](auto x){
        cout << *(x+1) << '\n';
        }, mt.variant_pointer);*/
//    std::cout << matrix.data.DATA4I[10] << std::endl;
    /*unsigned char* cp = new unsigned char [3 * sizeof(unsigned char)]{'a','b','c'};
    unsigned char* cp2 = new unsigned char [3 * sizeof(unsigned char)]{'a','b','a'};
    Matrix matrix1(Matrix::TYPE8U,1,3,cp,1,1);
//    Matrix matrix2(Matrix::TYPE8U,1,3,cp2,1,1);
    Matrix matrix2(matrix1);
    visit([](auto x){
        std::cout << *(x+1) << std::endl;
    },matrix2.variant_pointer);*/
    Matrix mt(Matrix::TYPE4I, 8, 8, ReadFromFile("afile.dat",Matrix::TYPE4I, 64));
    cout << mt << endl;

    /*Matrix mt1(Matrix::TYPE4I, 32, 8, ReadFromFile("afile.dat",Matrix::TYPE4I, 32*8));
    cout << mt1 << endl;

    Matrix mt2(Matrix::TYPE4I, 8, 32, ReadFromFile("afile.dat",Matrix::TYPE4I, 32*8));
    cout << mt2 << endl;*/

    mt.setElement(3,3,1,6.66);
    cout << mt << endl;
    float i;
    mt.getElement(3,3,1,i);
    cout << "i = " << i << endl;
    /*Matrix mtt(Matrix::TYPE64F);
    ROI(mt, mtt, 5, 5, 2, 2);
    cout << mtt << endl;*/
    return 0;
}