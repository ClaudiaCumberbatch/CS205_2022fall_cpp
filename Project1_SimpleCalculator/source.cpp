#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
using namespace std;

double *scientific_info(char *input)
{
    double *result = new double[3];
    long a1 = 0;   // a的整数位
    double a2 = 0; // a的小数位
    int k = 0; // 遍历字符数组的每个元素
    double dec = 0.1; // 计算小数
    bool flag = false; //标记是否来到小数点之后
    double a = 0;
    int digit = -1;
    while (input[k] != 0)
    {
        if (input[k] == '.')
        {
            flag = true; //来到小数点之后
            k++;
            continue;
        }
        else if (input[k] == 'e') //取exp
        {
            k++;
            int exp = 0;
            bool neg = false; //正数
            while (input[k] != 0)
            {
                if (input[k] == '-')
                {
                    neg = true;
                    k++;
                    continue;
                }
                else
                {
                    exp = exp * 10 + input[k] - '0';
                    k++;
                }
            }
            if (neg)
            {
                exp = -exp;
            }

            result[0] = a1;
            result[1] = a2;
            result[2] = exp;
            return result;
        }

        if (!flag)
        {
            digit++;
            a1 = a1 * 10 + input[k] - '0';
        }
        else
        {
            a2 = a2 + (input[k] - '0') * dec;
            dec = dec * 0.1;
        }
        k++;
    }
    a = a1 + a2;
    
    a1 = floor(a / pow(10, digit));
    a2 = (a - a1 * pow(10, digit)) * pow(0.1, digit);
    result[0] = a1;
    result[1] = a2;
    result[2] = digit;
    return result;
}

double scientific_num(double a1, double a2, double exp)
{
    double result = (a1 + a2) * pow(10, exp);
    return result;
}

bool legal(char *input)
{
    int k = 0;
    while (input[k] != 0)
    {
        if ((input[k] >= 30 && input[k] <= 57)                        
            || input[k] == '.' || input[k] == '-' 
            || input[k] == 'e' || input[k] == 'E') 
        {
            k++; //合法就往后挪
        }
        else
        {
            return false;
        }
    }
    return true;
}

double *get_num(char *input)
{
    double *output = new double[2];
    double a1 = 0;   // a的整数位
    double a2 = 0; // a的小数位
    int k = 0;
    double dec = 0.1;
    bool flag = false;
    
    while (input[k] != 0)
    {
        if (input[k] == '.')
        {
            flag = true; //来到小数点之后
            k++;
            continue;
        }
        else if (input[k] == 'e' || input[k] == 'E') //进入科学计数法
        {
            double *result = scientific_info(input);
            if (result[2] > 150)
            {
                output[0] = 1.0 / 0.0;
                output[1] = 1.0 / 0.0;
                return output;
            }
            else
            {
                double temp = scientific_num(result[0], result[1], result[2]);
                output[0] = floor(temp);
                output[1] = temp - output[0];
                return output;
            }
        }

        if (!flag)
        {
            a1 = a1 * 10 + input[k] - '0';
        }
        else
        {
            a2 = a2 + (input[k] - '0') * dec;
            dec = dec * 0.1;
        }
        k++;
    }
    output[0] = a1;
    output[1] = a2;
    return output;
}

string mul(double *inputs[])
{
    long a1 = (long)inputs[0][0] * (long)inputs[1][0];
    double a2 = inputs[0][0] * inputs[1][1] + inputs[0][1] * inputs[1][0] + inputs[0][1] * inputs[1][1];
    
    stringstream ss;
    ss << setprecision(15) << a2;
    string output;

    if (a2 > 0.0000001)
    {
        output = to_string(a1 + a2);
    }
    else if (a2 == 0)
    {
        output = to_string(a1);
    }
    else
    {
        output = to_string(a1) + " + " + ss.str();
    }
    return output;
}

string fuzzy_mul(char *inputs[])
{
    double *result1 = new double[3];
    double *result2 = new double[3];
    result1 = scientific_info(inputs[1]);
    result2 = scientific_info(inputs[2]);

    double a = (result1[0] + result1[1]) * (result2[0] + result2[1]);
    double exp = result1[2] + result2[2];
    if (a >=10)
    {
        a = a / 10;
        exp++;
    }
    string output = to_string(a) + "e" + to_string(exp);
    return output;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cerr << "Please input exactly two numbers!" << endl;
        return 0;
    }
    
    if (!(legal(argv[1]) && legal(argv[2])))
    {
        cerr << "The input cannot be interpret as numbers!" << endl;
        return 0;
    }
    else
    {
        double *a = get_num(argv[1]);
        double *b = get_num(argv[2]);

        if (isinf(a[0]) || isinf(b[0])) //超范围，进入fuzzy_mul
        {
            cout << "The input is so big that fuzzy multiplication is invoked!" << endl;
            cout << argv[1] << " * " << argv[2] << " = " << fuzzy_mul(argv) << endl;
        }
        else if (a[0] > 1e10 || b[0] > 1e10)
        {
            cout << argv[1] << " * " << argv[2] << " = " << fuzzy_mul(argv) << endl;
        }
        else
        {
            double *input[2] = {a, b};
            cout.precision(15);
            cout << argv[1] << " * " << argv[2] << " = " << mul(input) << endl;
        }
    }
    return 0;
}
