#include <iostream>
#include <cstring>
#include <stack>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <unordered_map>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;

const int N = 5000007;

const double PI = acos(-1);

int n, m;
int limit = 1;
int L; //二进制的位数
int R[N];

int num_of_exp = 0;
int left_brac_pos = -1;
int right_brac_pos = 0;
bool have_brac;
stack<char> opt; //操作符栈

struct Big_number
{
    int length;
    int scale;
    int digits[256];
    bool sign;
};

stack<Big_number> val; //操作数栈

struct variable
{
    string name;
    Big_number value;
};

struct Complex
{
    double x, y;
    Complex(double x = 0, double y = 0) : x(x), y(y) {}
} a[N], b[N];

Complex operator*(Complex J, Complex Q)
{
    //模长相乘，幅度相加
    return Complex(J.x * Q.x - J.y * Q.y, J.x * Q.y + J.y * Q.x);
}
Complex operator-(Complex J, Complex Q)
{
    return Complex(J.x - Q.x, J.y - Q.y);
}
Complex operator+(Complex J, Complex Q)
{
    return Complex(J.x + Q.x, J.y + Q.y);
}

//定义两个常量表示在数字中的状态
const int IN = 0;
const int OUT = 1;
variable var[256];
int flag = 0;
Big_number ans;

double get_num(string input);
int level(char theOpt);
bool change(string &from, string &to);
bool compute(string &theExp);
int have_var(string name);
string get_var_name(string exp);
bool is_numeric(char c);
Big_number input(string exp);
string output(Big_number big_number);
Big_number add(Big_number x, Big_number y);
Big_number miinus(Big_number x, Big_number y);
Big_number mul(Big_number x, Big_number y);
Big_number divide(Big_number x, Big_number y);
Big_number big_numbers[256];
int big_num_pos = 0; //用来访问big_numbers中的元素
void FFT(Complex *A, int type);
bool bigger(Big_number x, Big_number y);

int main()
{
    cout << "welcome to calculator 2.0" << endl;
    cout << "This is free software with ABSOLUTELY NO WARRANTY." << endl;
    ans = input("0");
    while (true)
    {
        //输入表达式
        string init_exp;
        cout << "------------------------" << endl;
        cout << "please input something you want to calculate:" << endl;
        cout << "or you can type \"q\" or \"quit\" to exit" << endl;
        cout << "------------------------" << endl;
        cin >> init_exp;
        num_of_exp++;
        have_brac = 0;
        right_brac_pos = init_exp.length();
        big_num_pos = 0;

        //sqrt
        if (init_exp.substr(0,4) == string("sqrt"))
        {
            int temp_length = init_exp.length() - 6;
            double num = get_num(init_exp.substr(5, temp_length));
            double left = 1;
            double right = num;
            double mid = 0.0;
            double result = 0.0;
            while (abs(result - num) > 0.0000000001)
            {
                mid = (left + right) / 2;
                result = mid * mid;
                if (result < num)
                {
                    left = mid;
                }
                else
                {
                    right = mid;
                }
            }
            ans = input(to_string(mid));
        }
        

        //判断是否要退出
        else if (init_exp == string("q") || init_exp == string("quit"))
        {
            cout << "------------------------" << endl;
            cout << "BYE" << endl;
            break;
        }

        else
        {
            for (int i = init_exp.length() - 1; i >= 0;)
            //括号，记录位置
            //变量，替换
            //等号，把括号里的东西计算掉
            //其他，接着往前挪
            {
                // cout << "this initial exp is " << init_exp << endl;
                // cout << "c " << i << " = " << init_exp[i] << endl;
                if (init_exp[i] == ')')
                {
                    right_brac_pos = i;
                    i--;
                    have_brac = 1;
                }
                else if (init_exp[i] == '(')
                {
                    left_brac_pos = i;
                    i--;
                    if (i < 0)
                    {
                        string infix = init_exp;
                        string postfix = "";
                        change(infix, postfix);
                        // cout << "postfix is " << postfix << endl;
                        compute(postfix);
                        ans = input(output(val.top()));
                    }
                    
                }
                else if (is_numeric(init_exp[i]))
                {
                    // cout << "init_exp is " << init_exp << endl;
                    // cout << "get numeric " << i << " " << init_exp[i] << endl;
                    if (i == 0) //只有计算，没有赋值
                    {
                        string infix = init_exp;
                        string postfix = "";
                        change(infix, postfix);
                        // cout << "postfix is " << postfix << endl;
                        compute(postfix);
                        ans = input(output(val.top()));
                        // cout << output(ans) << endl;
                    }
                    i--;
                }
                else if (init_exp[i] == '=' && init_exp[i - 1] != '=') //单个等号，此时等号后面已经没有变量
                // x=(y=1)赋值在括号内部，不知道左括号的位置
                // x=(1+1)赋值在括号外部，已经知道左括号位置，还没处理
                {
                    int length = (right_brac_pos - 1) - (i + 1) + 1; //表达式的长度
                    // cout << "right_bracket is " << right_brac_pos << endl;
                    // cout << "length is " << length << endl;
                    string infix = init_exp.substr(i + 1, length);
                    // cout << "infix is " << infix << endl;
                    string postfix = "";
                    change(infix, postfix);
                    // cout << "postfix is " << postfix << endl;
                    compute(postfix);

                    //找左括号，顺便获取变量的名字
                    string name = "";
                    left_brac_pos = -1;
                    for (int j = i - 1; j >= 0; j--)
                    {
                        if (init_exp[j] == '(')
                        {
                            left_brac_pos = j;
                            break;
                        }
                        else
                        {
                            name = init_exp[j] + name;
                        }
                    }
                    // cout << "name is " << name << endl;

                    //给变量赋值
                    int position = have_var(name);
                    if (position != -1) //原先有这个变量，更新他的值
                    {
                        var[position].value = input(output(val.top()));
                    }
                    else
                    {
                        var[flag].name = name;
                        var[flag].value = input(output(val.top()));
                        flag++;
                    }

                    //把整个括号里的东西替换掉
                    length = right_brac_pos - left_brac_pos + 1 - (1 - have_brac) * 2;
                    // cout << "length is " << length << endl;
                    Big_number temp = input(output(val.top()));
                    // cout << "left bracket position is " << left_brac_pos << endl;
                    // cout << "have_brac is " << have_brac << endl;
                    // cout << "length is " << length << endl;
                    // cout << "para is " << left_brac_pos + (1-have_brac) << endl;
                    init_exp.replace(left_brac_pos + (1 - have_brac), length, output(temp));
                    left_brac_pos = 0;
                    have_brac = 0;
                    right_brac_pos = init_exp.length();

                    //更新i
                    i = init_exp.length() - 1;
                    // cout << "init_exp is " << init_exp << endl;
                    // cout << "i = " << i << endl;
                }

                else //变量
                {
                    int head = i;
                    int tail = i;
                    for (int j = i; j >= 0; j--)
                    {
                        if (is_numeric(init_exp[j]))
                        {
                            head = j + 1;
                            break;
                        }
                    }
                    string name = get_var_name(init_exp.substr(head, tail - head + 1)); //开头在什么地方？
                    // cout << "init_exp.sub is " << init_exp.substr(head, tail - head + 1) << endl;
                    // cout << "name is (" << name << ")" << endl;
                    int position = have_var(name);
                    int length = i - head + 1;
                    // cout << "head is " << head << " and length is " << length << endl;
                    // cout << "position = " << position << endl;
                    if (position != -1) //原先有这个变量，获取他的值，替换字符串
                    {
                        init_exp.replace(head, length, output(var[position].value));
                        // cout << "exp is " << init_exp << endl;
                    }
                    else
                    {
                        init_exp.replace(head, length, "0");
                    }

                    //更新i
                    i = head - 1;
                    if (i < 0)
                    {
                        i = 0;
                    }

                    // cout << "i = " << i << endl;
                }
            }
        }

        //非法输入
        // else
        // {
        //     cout << "(standard_in) " << num_of_exp << " : parse error";
        // }

        //把变量全部打印一遍
        // for (int i = 0; i < 256; i++)
        // {
        //     if (var[i].name != "")
        //     {
        //         cout << var[i].name << " = " << output(var[i].value) << endl;
        //     }
        // }
        cout << output(ans) << endl;
    }
    return 0;
}

//字符串转化为数字
double get_num(string input)
{
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
    
    return a1 + a2;
}

//为每一个操作符返回一个数，数代表了优先级
int level(char theOpt)
{
    switch (theOpt)
    {
    case '(':
        return 0;
    case '+':
    case '-':
        return 1;
    case '*':
    case '/':
        return 2;
    case ')':
        return 3;
    }
    return -1;
}

//将中缀表达式转换成后缀表达式
bool change(string &from, string &to)
{
    int state = OUT;
    char c;

    for (int i = 0; i < from.length(); i++)
    {
        c = from[i];
        if (isdigit(c))
        { //是小数点前的数字
            to = to + c;
            state = IN; //状态改为在数字内
        }
        else
        {
            if (state == IN && c == '.')
            { //碰到小数点，仍旧是数字
                to = to + '.';
                continue;
            }
            if (state == IN && c != '.')
            { //不是数字
                to += ' ';
            }
            if (c == '=') 
                break;
            else if (c == '(')
                opt.push(c);
            else if (c == ')')
            {
                while (!opt.empty() && opt.top() != '(')
                { //括号匹配
                    to += opt.top();
                    to += ' ';
                    opt.pop();
                }
                opt.pop();
            }
            else
            {
                while (true)
                {
                    if (opt.empty() || opt.top() == '(')
                        opt.push(c);
                    else if (level(c) > level(opt.top()))
                        opt.push(c);
                    else
                    {
                        to += opt.top();
                        to += ' ';
                        opt.pop();
                        continue;
                    }
                    break;
                }
            }
            state = OUT; //状态为在数字外
        }
    }
    while (!opt.empty())
    {
        to += opt.top();
        to += ' ';
        opt.pop();
    }
    return true;
}

//计算后缀表达式
bool compute(string &theExp) //表达式中的数字转化成Big_number后压栈计算，加减法正常计算，乘法fft，除法用减法模拟
{
    int state = OUT; //初始状态为在数字外
    char c;
    bool single = true;
    int last_pos = -1; //上一个运算符的位置
    bool dot = false;

    int pos = 0;
    int scale = 0;
    int length = 0;
    for (int i = 0; i < theExp.length(); i++)
    {
        c = theExp[i];
        if (isdigit(c) || c == '.') //是数字
        {
            if (isdigit(c))
            {
                length++;
                big_numbers[big_num_pos].digits[pos] = c - '0';
                pos++;
                state = IN; //状态为在数字内
                if (dot == true)
                    scale++;
            }
            if (c == '.')
            {
                dot = true;
                continue;
            }
        }
        else //不是数字（或者小数点）
        {
            single = false;
            dot = false;
            if (state == IN)
            {
                big_numbers[big_num_pos].length = length;
                big_numbers[big_num_pos].scale = scale;
                big_numbers[big_num_pos].sign = true; //待定
                // cout << "ans = " << output(big_numbers[big_num_pos]) << endl;
                val.push(big_numbers[big_num_pos]);
                big_num_pos++;
                pos = 0;
                length = 0;
                scale = 0;
            }

            // 四则运算
            Big_number x, y;
            if (c != ' ')
            {
                x = val.top();
                // cout << "x = " << output(x) << endl;
                val.pop();
                y = val.top();
                // cout << "y = " << output(y) << endl;
                val.pop();

                switch (c)
                {
                case '+': //同号加法和异号减法 调用add
                    val.push(add(x, y));
                    single = false;
                    break;
                case '-':
                    val.push(miinus(y, x));
                    single = false;
                    break;
                case '*':
                    val.push(mul(x, y));
                    single = false;
                    break;
                case '/':
                    val.push(divide(y, x));
                    single = false;
                    break;
                default:
                    cout << "未知的错误!" << endl;
                }
            }

            state = OUT;
        }
    }
    //没有表达式的情况
    if (single)
    {
        if (state == IN)
        {
            big_numbers[big_num_pos].length = length;
            big_numbers[big_num_pos].scale = scale;
            big_numbers[big_num_pos].sign = true; //待定
            // cout << "ans = " << output(big_numbers[big_num_pos]) << endl;
            val.push(big_numbers[big_num_pos]);
            big_num_pos++;
            pos = 0;
            length = 0;
            scale = 0;
        }
    }

    return true;
}

Big_number add(Big_number x, Big_number y) //正数加法
{
    Big_number result;
    int x_before = x.length - x.scale;
    int y_before = y.length - y.scale;
    int length = max(x.scale, y.scale) + max(x_before, y_before) + 1; //最后需要更新成digits数组的大小
    int *cur_x = new int[length]();
    int *cur_y = new int[length]();
    int carry = 0;

    //整数部分补齐
    for (int i = max(x_before, y_before); i >= 0; i--) //填充前max(x_before, y_before)位
    {
        if (i <= max(x_before, y_before) - x_before)
        {
            cur_x[i] = 0;
        }
        else
        {
            int filled = max(x_before, y_before) - i;
            cur_x[i] = x.digits[x_before - filled - 1];
        }
        if (i <= max(x_before, y_before) - y_before)
        {
            cur_y[i] = 0;
        }
        else
        {
            int filled = max(x_before, y_before) - i;
            cur_y[i] = y.digits[y_before - filled - 1];
        }
        // cout << "cur_x " << i << " = " << cur_x[i] << endl;
        // cout << "cur_y " << i << " = " << cur_y[i] << endl;
    }

    //小数部分补齐
    for (int i = max(x.scale, y.scale) - 1; i >= 0; i--)
    {
        if (i + 1 > x.scale)
        {
            cur_x[max(x_before, y_before) + i + 1] = 0;
        }
        else
        {
            cur_x[max(x_before, y_before) + i + 1] = x.digits[x_before + i];
        }
        if (i + 1 > y.scale)
        {
            cur_y[max(x_before, y_before) + i + 1] = 0;
        }
        else
        {
            cur_y[max(x_before, y_before) + i + 1] = y.digits[y_before + i];
        }
        // cout << "cur_x " << max(x_before, y_before)+i+1 << " = " << cur_x[max(x_before, y_before)+i+1] << endl;
        // cout << "cur_y " << max(x_before, y_before)+i+1 << " = " << cur_y[max(x_before, y_before)+i+1] << endl;
    }

    for (int i = length - 1; i >= 0; i--)
    {
        result.digits[i] = (cur_x[i] + cur_y[i] + carry) % 10;
        carry = (cur_x[i] + cur_y[i] + carry) / 10;
    }

    if (result.digits[0] == 0) //处理最大位进位
    {
        for (int i = 0; i < length - 1; i++)
        {
            result.digits[i] = result.digits[i + 1];
        }
        length--;
    }

    result.length = length;
    result.scale = max(x.scale, y.scale);
    result.sign = x.sign; //存疑

    delete[] cur_x;
    delete[] cur_y;

    // cout << "result is " << output(result) << endl;
    // cout << "length is " << result.length << " and scale is " << result.scale << endl;

    return result;
}

Big_number miinus(Big_number x, Big_number y) // x-y且x>y，正数减法
{
    // cout << "x inside is " << output(x) << endl;
    // cout << "y inside is " << output(y) << endl;
    Big_number result;
    Big_number complement;
    int length = max(x.length, y.length);
    int com_pos = y.length - 1;
    int *cur_y = new int[length]();
    int x_before = x.length - x.scale;
    int y_before = y.length - y.scale;

    //小数部分补齐
    for (int i = max(x.scale, y.scale) - 1; i >= 0; i--)
    {
        if (i + 1 > y.scale)
        {
            cur_y[max(x_before, y_before) + i] = 0;
        }
        else
        {
            cur_y[max(x_before, y_before) + i] = y.digits[y_before + i];
        }
        // cout << "cur_y " << max(x_before, y_before)+i << " = " << cur_y[max(x_before, y_before)+i+1] << endl;
    }

    //整数部分补齐
    for (int i = max(x_before, y_before); i > 0; i--) //填充前max(x_before, y_before)位
    {
        if (i <= max(x_before, y_before) - y_before)
        {
            cur_y[i - 1] = 0;
        }
        else
        {
            int filled = max(x_before, y_before) - i;
            cur_y[i - 1] = y.digits[y_before - filled - 1];
        }
        // cout << "cur_y " << i << " = " << cur_y[i] << endl;
    }

    //计算complement
    for (int i = length - 1; i >= 0; i--)
    {
        if (com_pos >= 0)
        {
            complement.digits[i] = 9 - cur_y[i];
        }
        else
        {
            complement.digits[i] = 9;
        }
    }
    delete[] cur_y;
    complement.length = length;
    complement.scale = max(x.scale, y.scale);
    // complement.scale = 0;
    complement.sign = true;
    // cout << "complement is " << output(complement) << endl;
    result = add(x, complement);
    // cout << "big result is " << output(result) << endl;

    //去头
    for (int i = 0; i < length; i++)
    {
        result.digits[i] = result.digits[i + 1];
    }
    result.length = length;
    result.scale = 0;
    result = add(result, input("1"));
    result.scale = max(x.scale, y.scale);
    // cout << "length is " << result.length << endl;
    // cout << "small result is " << output(result) << endl;

    return result;
}

Big_number mul(Big_number x, Big_number y)
{
    m = x.length - 1;
    n = y.length - 1;
    Big_number result;
    int *temp = new int[m + n + 1]();
    for (int i = 0; i < x.length; i++)
    {
        a[i].x = x.digits[i];
    }
    for (int i = 0; i < y.length; i++)
    {
        b[i].x = y.digits[i];
    }
    while (limit <= n + m)
        limit <<= 1, L++;

    for (int i = 0; i < limit; i++)
    {
        R[i] = (R[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    FFT(a, 1); // FFT 把a的系数表示转化为点值表示
    FFT(b, 1); // FFT 把b的系数表示转化为点值表示

    //计算两个系数表示法的多项式相乘后的点值表示
    for (int i = 0; i <= limit; ++i)
        a[i] = a[i] * b[i];

    FFT(a, -1);

    // temp暂存多项式系数
    for (int i = 0; i <= m + n; i++)
    {
        // cout << "a" << i << " = " << a[i].x << endl; //错误在这之前
        temp[i] = (int)(a[i].x + 0.5);
        // cout << "temp " << i << " is " << temp[i] << endl;
    }

    //计算结果值
    int carry = 0;
    int temp_pos = m + n;
    for (int i = m + n + 2; i >= 0; i--)
    {
        if (temp_pos >= 0)
        {
            result.digits[i] = (carry + temp[temp_pos]) % 10;
            carry = (carry + temp[temp_pos]) / 10;
            temp_pos--;
        }
        else
        {
            result.digits[i] = carry % 10;
            carry = carry / 10;
        }
    }

    //去头
    int count = 0;
    while (result.digits[0] == 0)
    {
        count++;
        for (int i = 0; i < m + n + 3; i++)
        {
            result.digits[i] = result.digits[i + 1];
            // cout << "digit " << i << " = " << result.digits[i] << endl;
        }
    }

    result.scale = x.scale + y.scale;
    result.length = x.length + y.length + 1 - count;
    result.sign = !(x.sign ^ y.sign);
    // cout << "scale is " << result.scale << endl;
    // cout << "length is " << result.length << endl;
    // cout << "sign is " << result.sign << endl;
    // cout << "result inside is " << output(result) << endl;

    return result;
}

void FFT(Complex *A, int type)
{
    for (int i = 0; i < limit; ++i)
        if (i < R[i])
            swap(A[i], A[R[i]]);
    // i小于R[i]时才交换，防止同一个元素交换两次，回到它原来的位置。

    //从底层往上合并
    for (int mid = 1; mid < limit; mid <<= 1)
    {
        //待合并区间长度的一半，最开始是两个长度为1的序列合并,mid = 1;
        Complex wn(cos(PI / mid), type * sin(PI / mid)); //单位根w_n^1;

        for (int len = mid << 1, pos = 0; pos < limit; pos += len)
        {
            // len是区间的长度，pos是当前的位置,也就是合并到了哪一位
            Complex w(1, 0); //幂,一直乘，得到平方，三次方...

            for (int k = 0; k < mid; ++k, w = w * wn)
            {
                //只扫左半部分，蝴蝶变换得到右半部分的答案,w 为 w_n^k
                Complex x = A[pos + k];           //左半部分
                Complex y = w * A[pos + mid + k]; //右半部分
                A[pos + k] = x + y;               //左边加
                A[pos + mid + k] = x - y;         //右边减
            }
        }
    }
    if (type == 1)
        return;
    for (int i = 0; i <= limit; ++i)
        a[i].x /= limit;
    //最后要除以limit也就是补成了2的整数幂的那个N，将点值转换为系数
    //（前面推过了点值与系数之间相除是N）
}

Big_number divide(Big_number x, Big_number y)
{
    Big_number result;
    double x_num = 0;
    double y_num = 0;
    for (int i = 0; i < x.length; i++)
    {
        x_num = x_num * 10 + x.digits[i];
    }
    x_num = x_num / pow(10, x.scale);
    for (int i = 0; i < y.length; i++)
    {
        y_num = y_num * 10 + y.digits[i];
    }
    y_num = y_num / pow(10, y.scale);

    double temp = x_num / y_num;
    result = input(to_string(temp));

    return result;
}

// x是否大于等于y
bool bigger(Big_number x, Big_number y)
{
    if (x.sign == true && y.sign == false) // x正y负
    {
        return true;
    }
    else if (x.sign == false && y.sign == true) // x负y正
    {
        return false;
    }
    else if (x.sign == false) //都是正数
    {
        if (x.length - x.scale > y.length - y.scale)
        {
            return true;
        }
        else if (x.length - x.scale > y.length - y.scale)
        {
            return false;
        }
        else //整数部分长度相等
        {
            for (int i = x.length - x.scale - 1; i >= 0; i--)
            {
                if (x.digits[i] > y.digits[i])
                {
                    return true;
                }
                else if (x.digits[i] < y.digits[i])
                {
                    return false;
                }
            }
            for (int i = x.length - x.scale; i < min(x.length, y.length); i++)
            {
                if (x.digits[i] > y.digits[i])
                {
                    return true;
                }
                else if (x.digits[i] < y.digits[i])
                {
                    return false;
                }
            }
            if (x.length > y.length)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        for (int i = x.length - x.scale - 1; i >= 0; i--)
        {
            if (x.digits[i] > y.digits[i])
            {
                return false;
            }
            else if (x.digits[i] < y.digits[i])
            {
                return true;
            }
        }
        for (int i = x.length - x.scale; i < min(x.length, y.length); i++)
        {
            if (x.digits[i] > y.digits[i])
            {
                return false;
            }
            else if (x.digits[i] < y.digits[i])
            {
                return true;
            }
        }
        if (x.length > y.length)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    return true;
}

int have_var(string name)
{
    for (int i = 0; i < 256; i++) //这个遍历也太不优雅了吧
    {
        if (var[i].name == name)
        {
            return i;
        }
    }
    return -1;
}

string get_var_name(string exp)
{
    for (int i = 0; i < exp.size(); i++)
    {
        if (exp[i] == '=' || exp[i] == '+' || exp[i] == '-' || exp[i] == '*' || exp[i] == '/')
        {
            string output = exp.substr(0, i);
            return output;
        }
    }
    return exp;
}

bool is_numeric(char c)
{
    if (c == '+' || c == '-' || c == '*' || c == '/' || c == '.' || c == '(' || c == ')' || (c >= 48 && c <= 57))
    {
        return true;
    }
    else
    {
        return false;
    }
}

Big_number input(string exp) //并且去掉小数点后末尾的0
{
    Big_number result;
    int length = exp.length();
    int scale = 0;
    int l = 0;
    result.digits[exp.length() - 1] = 0;
    if (exp[0] == '-')
    {
        result.sign = 0;
    }
    else
    {
        result.sign = 1;
    }

    for (int i = 0; i < exp.length(); i++)
    {
        if (exp[i] == '.')
        {
            scale = exp.length() - i - 1;
            length--;
        }
        else
        {
            result.digits[l] = exp[i] - '0';
            l++;
        }
    }

    if (scale != 0)
    {
        for (int i = l - 1; i > 0; i--)
        {
            if (result.digits[i] == 0)
            {
                length--;
                scale--;
            }
            else
            {
                break;
            }
        }
    }

    result.length = length;
    result.scale = scale;

    return result;
}

string output(Big_number big_number)
{
    string output;
    int head_pos = 0;
    if (big_number.sign) //正数
    {
        output = "";
    }
    else
    {
        output = "-";
    }

    for (int i = 0; i < big_number.length; i++) // 去掉开头的0
    {
        if (i == big_number.length - big_number.scale)
        {
            output = output + ".";
        }
        char c = big_number.digits[i] + '0';
        output = output + c;
    }
    return output;
}