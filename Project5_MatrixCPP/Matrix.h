//
// Created by 周思呈 on 2022/12/11.
//

#ifndef PROJECT_5_MATRIXCPP_MATRIX_H
#define PROJECT_5_MATRIXCPP_MATRIX_H

#include <cstddef>
#include <memory>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <fstream>

/*ERROR TABLE*/
#define ERROR_TABLE(eid) if (eid == 1) std::cerr << "Illegal type in " << __FILE__ << ": Line " << __LINE__ << " in function " << __FUNCTION__ << std::endl;\
                        else if (eid == 2) std::cerr << "Illegal size in " << __FILE__ << ": Line " << __LINE__ << " in function " << __FUNCTION__ << std::endl;\
                        else if (eid == 3) std::cerr << "Index out of range in " << __FILE__ << ": Line " << __LINE__ << " in function " << __FUNCTION__ << std::endl;\
                        else if (eid == 4) std::cerr << "Empty matrix in " << __FILE__ << ": Line " << __LINE__ << " in function " << __FUNCTION__ << std::endl;\
                        else if (eid == 5) std::cerr << "Null pointer in " << __FILE__ << ": Line " << __LINE__ << " in function " << __FUNCTION__ << std::endl;\
                        else std::cerr << "Unknown error in " << __FILE__ << ": Line " << __LINE__ << " in function " << __FUNCTION__ << std::endl;
#define PRECISION 5
#define WIDTH 5

/*class ptr_manager{ //构造函数输入一个数组指针，创建一个shared_ptr成员
public:
    union pointers
    {
        std::shared_ptr<int*> dp;
    };
    std::shared_ptr<int*> dp;
    ptr_manager(){
        std::cout << "default constructor is invoked" << std::endl;
    }
    ptr_manager(int* p){
        std::cout << "ptr_manager(int* p) is invoked____" << std::endl;
        if (p == NULL) std::cout << "null pointer" << std::endl;
        dp = std::make_shared<int*>(p);
    }
};*/

class Matrix {
public:
    enum DataType {
        TYPE8U,
        TYPE8S,
        TYPE4I,
        TYPE32F,
        TYPE64F
    };
    DataType TYPE;
    size_t SIZE;
    /*union DATA
    {
        void * DATANULL;
        unsigned char * DATA8U;
        short *DATA8S;
        int *DATA4I;
        float *DATA32F;
        double *DATA64F;
    }data;*/
    std::variant<unsigned char *, short *, int *, float *, double *> variant_pointer;
    size_t ROW;
    size_t COL;
    size_t STEP;
    int CHANNEL;
    size_t *ref_cnt;

    Matrix(DataType type = TYPE8U, size_t row = 1, size_t col = 1, void *dp = NULL, int channel = 1, size_t step = 0);

    Matrix(const Matrix &m);

    bool release();

    ~Matrix();

    Matrix &operator=(const Matrix &m);

    bool operator==(const Matrix &m);

    bool operator!=(const Matrix &m);

    Matrix operator+(const Matrix &m) const;

    template <typename T>
    Matrix operator+(T c) const
    {
        try{
            if (ROW == 0 || COL == 0 || CHANNEL == 0 || STEP == 0)
                throw 4;

            size_t size = ROW * COL * CHANNEL;
            Matrix sum(TYPE, ROW, COL, NULL, CHANNEL, STEP);
            visit([&size, &sum, &c, this](auto x){
                visit([&size, &sum, &x, &c, this](auto y){
                    if (!x || !y) {
                        sum.release();
                        throw 4;
                    }
                    for (size_t i = 0; i < ROW; i++){
                        size_t temp = i * STEP;
                        for (size_t j = 0; j < c; j++){
                            *(y + temp + j) = *(x + temp + j) + c;
                        }
                    }
                },sum.variant_pointer);
            }, variant_pointer);
            return sum;
        }catch(int eid){
            ERROR_TABLE(eid)
        }
        return *this;
    }

    template <typename T>
    friend Matrix operator+(T c, const Matrix & m)
    {
        return m+c;
    }

    Matrix operator-(const Matrix &m) const;

    template <typename T>
    Matrix operator-(T c) const
    {
        return *this + (-c);
    }

    template <typename T>
    friend Matrix operator-(T c, const Matrix & m)
    {
        return m-c;
    }

    Matrix operator*(const Matrix & m) const;

    template <typename T>
    Matrix operator*(T c) const
    {
        try{
            if (ROW == 0 || COL == 0 || CHANNEL == 0 || STEP == 0)
                throw 4;

            size_t size = ROW * COL * CHANNEL;
            Matrix sum(TYPE, ROW, COL, NULL, CHANNEL, STEP);
            visit([&size, &sum, &c, this](auto x){
                visit([&size, &sum, &x, &c, this](auto y){
                    if (!x || !y) {
                        sum.release();
                        throw 4;
                    }
                    for (size_t i = 0; i < ROW; i++){
                        size_t temp = i * STEP;
                        for (size_t j = 0; j < c; j++){
                            *(y + temp + j) = *(x + temp + j) * c;
                        }
                    }
                },sum.variant_pointer);
            }, variant_pointer);
            return sum;
        }catch(int eid){
            ERROR_TABLE(eid)
        }
        return *this;
    }

    template <typename T>
    friend Matrix operator*(T c, const Matrix & m)
    {
        return m * c;
    }

    template <typename T>
    Matrix operator/(T c) const
    {
        return (*this) * (1/double(c));
    }

    friend std::ostream & operator<<(std::ostream & os, const Matrix & m);

    friend std::ostream & operator<<(std::ostream & os, const Matrix::DataType &tp);

    template <typename T>
    bool setElement(size_t r, size_t c, int channel, T value) {
        try {
            if (r >= ROW || c >= COL || channel > CHANNEL)
                throw 3;
            size_t index = r * CHANNEL * STEP + c * channel;
            visit([&index, &value](auto x){
                if (!x) throw 4;
                *(x + index) = value;
            }, variant_pointer);
            return true;
        } catch (int eid) {
            ERROR_TABLE(eid)
        }
        return false;
    }

    template <typename T>
    bool getElement(size_t r, size_t c, int channel, T &value) {
        try {
            if (r >= ROW || c >= COL || channel > CHANNEL)
                throw 3;
            size_t index = r * CHANNEL * STEP + c * channel;
            visit([&index, &value](auto x){
                if (!x) throw 4;
                value = *(x + index);
            }, variant_pointer);
            return true;
        } catch (int eid) {
            ERROR_TABLE(eid)
        }
        return false;
    }
};

bool ROI(const Matrix& original, Matrix &result, size_t indexR = 0, size_t indexC = 0, size_t rows = 0, size_t cols = 0);
Matrix deepCopy(const Matrix& original, size_t indexR = 0, size_t indexC = 0, size_t rows = 0, size_t cols = 0);
void* ReadFromFile(std::string filename, Matrix::DataType type, size_t size);

#endif //PROJECT_5_MATRIXCPP_MATRIX_H
