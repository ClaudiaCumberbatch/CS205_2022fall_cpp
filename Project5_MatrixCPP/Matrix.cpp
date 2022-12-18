//
// Created by 周思呈 on 2022/12/11.
//

#include "Matrix.h"

Matrix::Matrix(DataType type, size_t row, size_t col, void* dp, int channel, size_t step) {
    try {
        //一通参数检查
        if (type > 4 || type < 0)
            throw 1;
        if (channel < 0)
            throw 2;
        if (row == 0 || col == 0)
            throw 4;

        TYPE = type;
        ROW = row;
        COL = col;
        CHANNEL = channel;
        if (step) STEP = step;
        else STEP = COL;
        ref_cnt = new size_t[sizeof(size_t)]{1};
        if (dp) { //传入非空指针
            switch (type) {
                case TYPE8U:
//                        data.DATA8U = static_cast<unsigned char *>(dp);
                    variant_pointer = (unsigned char *) dp;
                    SIZE = sizeof(unsigned char);
                    break;
                case TYPE8S:
//                        data.DATA8S = static_cast<short *>(dp);
                    variant_pointer = (short *) dp;
                    SIZE = sizeof(short);
                    break;
                case TYPE4I:
//                        data.DATA4I = static_cast<int *>(dp);
                    variant_pointer = (int *) dp;
                    SIZE = sizeof(int);
                    break;
                case TYPE32F:
//                        data.DATA32F = static_cast<float *>(dp);
                    variant_pointer = (float *) dp;
                    SIZE = sizeof(float);
                    break;
                case TYPE64F:
//                        data.DATA64F = static_cast<double *>(dp);
                    variant_pointer = (double *) dp;
                    SIZE = sizeof(double);
                    break;
                default:
                    throw 666;
            }
        } else { //指针为空，需要申请内存
            switch (type) {
                case TYPE8U:
//                        data.DATA8U = new unsigned char[row * col * channel * sizeof(unsigned char)];
                    dp = new unsigned char[row * col * channel * sizeof(unsigned char)];
                    variant_pointer = (unsigned char *) dp;
                    SIZE = sizeof(unsigned char);
                    break;
                case TYPE8S:
//                        data.DATA8S = new short [row * col * channel * sizeof(short)];
                    dp = new short[row * col * channel * sizeof(short)];
                    variant_pointer = (short *) dp;
                    SIZE = sizeof(short);
                    break;
                case TYPE4I:
//                        data.DATA4I = new int [row * col * channel * sizeof(int)];
                    dp = new int[row * col * channel * sizeof(int)];
                    variant_pointer = (int *) dp;
                    SIZE = sizeof(int);
                    break;
                case TYPE32F:
//                        data.DATA32F = new float[row * col * channel * sizeof(float)];
                    dp = new float[row * col * channel * sizeof(float)];
                    variant_pointer = (float *) dp;
                    SIZE = sizeof(float);
                    break;
                case TYPE64F:
//                        data.DATA64F = new double[row * col * channel * sizeof(double)];
                    dp = new double[row * col * channel * sizeof(double)];
                    variant_pointer = (double *) dp;
                    SIZE = sizeof(double);
                    break;
                default:
                    std::cerr << "memory allocate fail" << std::endl;
                    throw 4;
            }
        }
    } catch (std::bad_alloc &ba) {
        std::cout << "bad_alloc exception!" << std::endl;
        std::cout << ba.what() << std::endl;
    } catch (int eid) {
        delete(ref_cnt);
        ERROR_TABLE(eid)
    }
}

Matrix::Matrix(const Matrix &m){
    try {
        //不存在空引用所以不用检查m是否为NULL
        if (m.ROW == 0 || m.COL == 0 || m.CHANNEL == 0 || m.STEP == 0) throw 4;
        if (!m.ref_cnt) throw 5;
        TYPE = m.TYPE;
        ROW = m.ROW;
        COL = m.COL;
        SIZE = m.SIZE;
        m.ref_cnt[0]++;
        ref_cnt = m.ref_cnt;
        visit([this](auto x){
            if (!x) throw 5;
            else variant_pointer = x;
        }, m.variant_pointer);
        variant_pointer = m.variant_pointer;
        CHANNEL = m.CHANNEL;
        STEP = m.STEP;
    }catch (int eid){
        ERROR_TABLE(eid)
    }
}

bool Matrix::release(){
    visit([this](auto x){
        if (!x) {
            return false;
        }else {
            delete(x);
            delete(ref_cnt);
            return true;
        }
    }, variant_pointer);
    return true;
}

Matrix::~Matrix(){
    try{
        if (!ref_cnt) throw 5;
        (*ref_cnt)--;
        if (!(*ref_cnt)){ //ref_cnt==0释放内存
            if (!this->release()) throw 5;
        }
    }catch(int eid){
        ERROR_TABLE(eid)
    }
}

Matrix &Matrix::operator=(const Matrix &m){
    try {
        if (m.ROW == 0 || m.COL == 0 || m.CHANNEL == 0 || m.STEP == 0) throw 4;
        if (!m.ref_cnt) throw 4;
        TYPE = m.TYPE;
        ROW = m.ROW;
        COL = m.COL;
        SIZE = m.SIZE;
        m.ref_cnt[0]++;
        ref_cnt = m.ref_cnt;
        visit([this](auto x){
            if (!x) throw 4;
            else variant_pointer = x;
        }, m.variant_pointer);
        CHANNEL = m.CHANNEL;
        STEP = m.STEP;
        return *this;
    }catch (int eid){
        ERROR_TABLE(eid)
    }
    return *this;
}

bool Matrix::operator==(const Matrix &m){
    if (TYPE != m.TYPE
        || ROW != m.ROW
        || COL != m.COL
        || CHANNEL != m.CHANNEL
        || STEP != m.STEP
        || SIZE != m.SIZE) {
        return false;
    }
    try{
        size_t size = ROW * COL * CHANNEL;
        switch (TYPE) {
            case TYPE8U:
            case TYPE8S:
            case TYPE4I:
                visit([&m](auto x){
                    visit([&x, &m](auto y){
                        if (!x || !y) throw 5;
                        for (size_t i = 0; i < m.ROW; i++){
                            size_t temp = i * m.STEP;
                            if (!memcmp(y+temp, x+temp, m.COL * m.CHANNEL * m.SIZE)) {
                                return false;
                            }else {
                                continue;
                            }
                        }
                    },m.variant_pointer);
                }, variant_pointer);
                break;
            case TYPE32F:
                visit([&m](auto x){
                    visit([&x, &m](auto y){
                        if (!x || !y) throw 5;
                        size_t c = m.COL * m.CHANNEL;
                        for (size_t i = 0; i < m.ROW; i++){
                            for (size_t j = 0; j < c; j++){
                                if (abs(*(x + i * m.STEP + j) - *(y + i * m.STEP + j)) > FLT_EPSILON) {
                                    return false;
                                }else{
                                    continue;
                                }
                            }
                        }
                        return true;
                    },m.variant_pointer);
                }, variant_pointer);
                break;
            case TYPE64F:
                visit([&m](auto x){
                    visit([&x, &m](auto y){
                        if (!x || !y) throw 5;
                        size_t c = m.COL * m.CHANNEL;
                        for (size_t i = 0; i < m.ROW; i++){
                            for (size_t j = 0; j < c; j++){
                                if (abs(*(x + i * m.STEP + j) - *(y + i * m.STEP + j)) > DBL_EPSILON) {
                                    return false;
                                }else{
                                    continue;
                                }
                            }
                        }
                        return true;
                    },m.variant_pointer);
                }, variant_pointer);
                break;
            default:
                throw 666;
        }
        return true;
    }catch (int eid){
        ERROR_TABLE(eid)
    }
    return false;
}

bool Matrix::operator!=(const Matrix &m){
    if (*this == m) return true;
    else return false;
}

Matrix Matrix::operator+(const Matrix & m) const
{
    try{
        if (ROW != m.ROW || COL != m.COL || CHANNEL != m.CHANNEL || STEP != m.STEP)
            throw 2;

        size_t c = COL * CHANNEL;
        DataType s_type = TYPE8U;
        if (TYPE <= m.TYPE){ //m的数据类型优先级更高
            s_type = m.TYPE;
        }else{
            s_type = TYPE;
        }
        Matrix sum(s_type, ROW, COL, NULL, CHANNEL, STEP);
        visit([&sum, &m, &c](auto x){
            visit([&sum, &x, &m, &c](auto y){
                visit([&x, &y, &m, &c, &sum](auto s){
                    if (!x || !y || !s) {
                        sum.release();
                        throw 5;
                    }
                    for (size_t i = 0; i < m.ROW; i++){
                        size_t temp = i * m.STEP;
                        for (size_t j = 0; j < c; j++){
                            *(s + temp + j) = *(x + temp + j) + *(y + temp + j);
                        }
                    }
                }, sum.variant_pointer);
            },m.variant_pointer);
        }, variant_pointer);
        return sum;
    }catch(int eid){
        ERROR_TABLE(eid)
    }
    return *this;
}

Matrix Matrix::operator-(const Matrix & m) const
{
    try{
        if (ROW != m.ROW || COL != m.COL || CHANNEL != m.CHANNEL || STEP != m.STEP)
            throw 2;

        size_t c = COL * CHANNEL;
        DataType s_type = TYPE8U;
        if (TYPE <= m.TYPE){ //m的数据类型优先级更高
            s_type = m.TYPE;
        }else{
            s_type = TYPE;
        }
        Matrix sum(s_type, ROW, COL, NULL, CHANNEL, STEP);
        visit([&sum, &m, &c](auto x){
            visit([&sum, &x, &m, &c](auto y){
                visit([&x, &y, &m, &c, &sum](auto s){
                    if (!x || !y || !s) {
                        sum.release();
                        throw 5;
                    }
                    for (size_t i = 0; i < m.ROW; i++){
                        size_t temp = i * m.STEP;
                        for (size_t j = 0; j < c; j++){
                            *(s + temp + j) = *(x + temp + j) - *(y + temp + j);
                        }
                    }
                }, sum.variant_pointer);
            },m.variant_pointer);
        }, variant_pointer);
        return sum;
    }catch(int eid){
        ERROR_TABLE(eid)
    }
    return *this;
}

Matrix Matrix::operator*(const Matrix & m) const
{
    try{
        if (ROW != m.COL || CHANNEL != m.CHANNEL)
            throw 2;
        size_t row = ROW;
        size_t col = m.COL;
        size_t mul = COL;
        DataType s_type = TYPE8U;
        if (TYPE <= m.TYPE){ //m的数据类型优先级更高
            s_type = m.TYPE;
        }else{
            s_type = TYPE;
        }
        Matrix sum(s_type, row, col, NULL, CHANNEL, m.STEP);
        visit([&row, &col, &mul, &sum, &m, this](auto x){
            visit([&row, &col, &mul, &sum, &x, &m, this](auto y){
                visit([&row, &col, &mul, &x, &y, &sum, &m, this](auto s){
                    if (!x || !y || !s) {
                        sum.release();
                        throw 4;
                    }
                    memset(s, 0, row * col * sum.SIZE);
#pragma omp parallel for
                    for (size_t i = 0; i < row; i++){
                        for (size_t k = 0; k < mul; k++){
                            for (size_t j = 0; j < col; j++){
//                                *(s + i * col + j) += *(x + i * mul + k) * *(y + k * col + j);
                                *(s + i * m.STEP + j) += *(x + i * STEP + k) * *(y + k * m.STEP + j);
                            }
                        }
                    }
                }, sum.variant_pointer);
            },m.variant_pointer);
        }, variant_pointer);
        return sum;
    }catch(int eid){
        ERROR_TABLE(eid)
    }
    return *this;
}

std::ostream & operator<<(std::ostream & os, const Matrix & m)
{
    try {
        size_t c = m.COL * m.CHANNEL;
        visit([&c, &os, &m](auto x) {
            if (!x) throw 4;
            if (c >= 16 && m.ROW >= 16) {
                os << "total " << m.ROW << " rows and " << c << " cols" << std::endl;
                for (int i = 0; i < 4; i++) { //前四行
                    for (int j = 0; j < 4; j++) { //前四列
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                           << *(x + i * m.STEP + j);
                    }
                    os << "......  ";
                    for (size_t j = c-4; j < c; j++) { //后四列
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                        << *(x + i * m.STEP + j);
                    }
                    os << std::endl;
                }
                os << "......" << std::endl;
                for (int i = m.ROW-4; i < m.ROW; i++) { //后四行
                    for (int j = 0; j < 4; j++) { //前四列
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                           << *(x + i * m.STEP + j);
                    }
                    os << "......  ";
                    for (size_t j = c-4; j < c; j++) { //后四列
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                           << *(x + i * m.STEP + j);
                    }
                    os << std::endl;
                }
            }else if (c >= 16){
                os << "total " << c << " cols" << std::endl;
                for (size_t i = 0; i < m.ROW; i++){
                    for (int j = 0; j < 4; j++) { //前四列
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                           << *(x + i * m.STEP + j);
                    }
                    os << "......  ";
                    for (size_t j = c-4; j < c; j++) { //后四列
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                           << *(x + i * m.STEP + j);
                    }
                    os << std::endl;
                }
            }else if (m.ROW >= 16){
                os << "total " << m.ROW << " rows" << std::endl;
                for (size_t i = 0; i < 4; i++){
                    for (size_t j = 0; j < c; j++){
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                           << *(x + i * m.STEP + j);
                    }
                    os << std::endl;
                }
                os << "......" << std::endl;
                for (size_t i = m.ROW-4; i < m.ROW; i++){
                    for (size_t j = 0; j < c; j++){
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                           << *(x + i * m.STEP + j);
                    }
                    os << std::endl;
                }
            }else{
                for (size_t i = 0; i < m.ROW; i++){
                    for (size_t j = 0; j < c; j++){
                        os << std::setprecision(PRECISION) << std::setw(WIDTH) << std::left
                           << *(x + i * m.STEP + j);
                    }
                    os << std::endl;
                }
            }
        }, m.variant_pointer);
        return os;
    }catch (int eid){
        ERROR_TABLE(eid)
    }
    return os;
}

std::ostream & operator<<(std::ostream & os, const Matrix::DataType &tp)
{
    switch (tp){
        case Matrix::TYPE8U:
            os << "TYPE8U";
            break;
        case Matrix::TYPE64F:
            os << "TYPE64F";
            break;
        case Matrix::TYPE8S:
            os << "TYPE8S";
            break;
        case Matrix::TYPE32F:
            os << "TYPE32F";
            break;
        case Matrix::TYPE4I:
            os << "TYPE4I";
            break;
        default:
            break;
    }
    return os;
}

bool ROI(const Matrix& original, Matrix &result, size_t indexR, size_t indexC, size_t rows, size_t cols){
    try{
        if (rows == 0) rows = original.ROW;
        if (cols == 0) cols = original.COL;
        if (original.ROW < indexR || indexR + rows >= original.ROW
        || original.COL < indexC || indexC + cols >= original.COL)
            throw 3;
        result.STEP = original.STEP;
        result.SIZE = original.SIZE;
        result.TYPE = original.TYPE;
        result.CHANNEL = original.CHANNEL;
        result.COL = cols;
        result.ROW = rows;
        size_t bias = original.COL * original.CHANNEL * indexR + indexC;
        visit([bias, &result](auto x){
            if (!x) throw 4;
            result.variant_pointer = x+bias;
        }, original.variant_pointer);
        if (!original.ref_cnt) throw 5;
        else (*original.ref_cnt)++;
        if (result.ref_cnt) delete(result.ref_cnt);
        result.ref_cnt = original.ref_cnt;
        return true;
    }catch(int eid){
        ERROR_TABLE(eid)
    }
    return false;
}

void* ReadFromFile(std::string filename, Matrix::DataType type, size_t size){
    try {
        char data[100]; //字符串
        int index = 0;
        // 以读模式打开文件
        std::ifstream infile;
        infile.open(filename);
//        std::cout << "Reading from the file" << std::endl;
        switch (type) {
            case Matrix::TYPE8U: {
                unsigned char *dp = new unsigned char[size * sizeof(unsigned char)];
                while (infile >> data && size > index) {
                    *(dp + index) = (unsigned char) (std::stoi(data));
                    index++;
                }
                infile.close();
                return dp;
            }
            case Matrix::TYPE8S: {
                short *dp = new short[size * sizeof(short)];
                while (infile >> data && size > index) {
                    *(dp + index) = short(std::stoi(data));
                    index++;
                }
                infile.close();
                return dp;
            }
            case Matrix::TYPE4I: {
                int *dp = new int[size * sizeof(int)];
                while (infile >> data && size > index) {
                    *(dp + index) = std::stoi(data);
                    index++;
                }
                infile.close();
                return dp;
            }
            case Matrix::TYPE32F:{
                float *dp = new float[size * sizeof(float)];
                while (infile >> data && size > index) {
                    *(dp + index) = std::stof(data);
                    index++;
                }
                infile.close();
                return dp;
            }
            case Matrix::TYPE64F:{
                double *dp = new double[size * sizeof(double)];
                while (infile >> data && size > index) {
                    *(dp + index) = std::stod(data);
                    index++;
                }
                infile.close();
                return dp;
            }
            default:
                throw 666;
        }
    }catch (std::bad_alloc & ba) {
        std::cout << "bad_alloc exception!" << std::endl;
        std::cout << ba.what() << std::endl;
    }catch (int eid){
        ERROR_TABLE(eid)
    }
}

Matrix deepCopy(const Matrix& original, size_t indexR, size_t indexC, size_t rows, size_t cols){
    try{
        if (rows == 0) rows = original.ROW;
        if (cols == 0) cols = original.COL;
        if (original.ROW < indexR || indexR + rows >= original.ROW
            || original.COL < indexC || indexC + cols >= original.COL)
            throw 3;
        Matrix result(original.TYPE, rows, cols, NULL, original.CHANNEL, original.STEP);
        visit([&result](auto x){
            visit([&x, &result](auto y){
                if (!x || !y) throw 5;
                for (size_t i = 0; i < result.ROW; i++){
                    size_t temp = i * result.STEP;
                    memcpy(y+temp, x+temp, result.COL * result.CHANNEL * result.SIZE);
                }
            },result.variant_pointer);
        }, original.variant_pointer);
        return result;
    }catch(int eid){
        ERROR_TABLE(eid)
    }
    return original;
}
