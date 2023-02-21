#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <unordered_map>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;

const int N = 5000007;

const double PI = acos(-1);

int n, m;
int res, ans[N];
int limit = 1;
int L; //二进制的位数
int R[N];

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

int main()
{
    scanf("%d%d", &n, &m);
    //读入多项式的每一项，保存在复数的实部

    for (int i = 0; i <= n; ++i) //输入第一个多项式
    {
        double x;
        scanf("%lf", &x);
        a[i].x = x; // complex类型变量.real(x)意味着将实数部赋为x，real()返回实数部值
    }
    for (int i = 0; i <= m; ++i) //输入第二个多项式
    {
        double x;
        scanf("%lf", &x);
        b[i].x = x;
    }
    while (limit <= n + m)
        limit <<= 1, L++;
    //也可以写成：limit = 1 << int(log2(n + m) + 1);
    cout << R[0] << " " << R[1] << endl;
    // 补成2的整次幂，也就是N
    for (int i = 0; i < limit; ++i)
    {
        R[i] = (R[i >> 1] >> 1) | ((i & 1) << (L - 1));
        cout << "R[i] = " << R[i] << endl;
    }
    FFT(a, 1); // FFT 把a的系数表示转化为点值表示
    FFT(b, 1); // FFT 把b的系数表示转化为点值表示
    //计算两个系数表示法的多项式相乘后的点值表示
    
    for (int i = 0; i <= limit; ++i)
        a[i] = a[i] * b[i];
    //对应项相乘，O(n)得到点值表示的多项式的解C，利用逆变换完成插值得到答案C的点值表示
    FFT(a, -1);

    for (int i = 0; i <= n + m; ++i)
        //这里的 x 和 y 是 double 的 hhh
        printf("%d ", (int)(a[i].x + 0.5)); //注意要+0.5，否则精度会有问题
}
