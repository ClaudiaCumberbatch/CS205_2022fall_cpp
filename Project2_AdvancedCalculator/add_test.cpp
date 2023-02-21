#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
using namespace std;

typedef vector<int> vint;

vint A, B, a, b; // 大写为整数部分，小写为小数部分

vint add(vint a, vint b) // 整数相加模板
{
    vint c;
    int t = 0;
    for (int i = 0; i < a.size() || i < b.size(); i++)
    {
        if (i < a.size())
            t += a[i];
        if (i < b.size())
            t += b[i];
        c.push_back(t % 10);
        t /= 10;
    }
    if (t)
        c.push_back(t);
    return c;
}

int main()
{
    string n, m;
    cin >> n >> m;
    bool flag = true;                       // 表示是否应该开始读取整数部分
    for (int i = n.size() - 1; i >= 0; i--) // 读第一个
    {
        if (n[i] == '.')
        {
            flag = false;
            continue;
        }
        if (flag)
            a.push_back(n[i] - '0');
        else
            A.push_back(n[i] - '0');
    }
    flag = true;
    for (int i = m.size() - 1; i >= 0; i--) // 读第二个
    {
        if (m[i] == '.')
        {
            flag = false;
            continue;
        }
        if (flag)
            b.push_back(m[i] - '0');
        else
            B.push_back(m[i] - '0');
    }

    while (a.size() != b.size()) // 小数部分补齐
    {
        if (a.size() > b.size())
            b.insert(b.begin(), 0);
        if (a.size() < b.size())
            a.insert(a.begin(), 0);
    }

    vint C = add(A, B), c = add(a, b);

    if (c.size() > a.size()) // 说明小数部分向整数进位
    {
        vint temp;
        temp.push_back(c[c.size() - 1]);
        C = add(C, temp);
        c.erase(c.end() - 1); // 删去进位
    }

    while (c.size() > 1 && c.front() == 0)
        c.erase(c.begin()); // 删除小数的尾0

    for (int i = C.size() - 1; i >= 0; i--)
        cout << C[i];
    cout << ".";
    for (int i = c.size() - 1; i >= 0; i--)
        cout << c[i];
    return 0;
}
