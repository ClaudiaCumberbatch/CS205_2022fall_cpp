#include <iostream>
#include <cstring>
#include <cmath>
#include <stack>
#include <iomanip>
#include <cstdlib>
#include <complex>
#include <complex.h>

using namespace std;

stack<char> opt;   //操作符栈
stack<double> val; //操作数栈

//求阶乘所用数组
int LargeNumberFactorial[1000] = {0};

//定义两个常量表示在数字中的状态
const int IN = 0;
const int OUT = 1;

char opt_set[10] = "+-*/()=.";

int level(char theOpt);								//操作符的优先级
bool in_set(char theChar);							//判断输入的操作符是否合法，仅包括()+-*/
bool del_space(string &theString);					//用于去除空格并检查是否有非法字符
string to_string(int theInt);						//把数字转换成字符串
bool change(string &from, string &to);				//将中缀表达式转换成后缀表达式
bool compute(string &theExp);						//计算后缀表达式
int Judge_rectangle(string &theStr);				//判断是否为三角函数或反三角函数并返回对应的编号
double cal_rect_number(string &theStr, int theInt); //求三角函数里包含的数字并求出对应的结果
bool Judge_factorial(string &theStr);				//判断是否为求阶乘
int cal_factorial(int theInt);						//利用数组求阶乘
int cal_factorial_number(string &theStr);			// 求阶乘所求的数字

int main()
{
	cout << "欢迎使用计算器" << endl;
	cout << "你可以使用这个计算器做一些简单的运算" << endl;
	cout << "如果你想退出，输入:exit即可..." << endl;
	while (true)
	{

		//输入表达式
		string init_exp;
		cout << "---请输入表达式:";
		cin >> init_exp;

		//判断是否要退出
		if (init_exp == string("exit"))
		{
			cout << "感谢您的使用..." << endl;
			break;
		}

		//判断是否三角函数和反三角函数并求值
		else if (Judge_rectangle(init_exp) != 0)
		{
			double res = 0;
			res = cal_rect_number(init_exp, Judge_rectangle(init_exp));
			if (res)
				cout << "结果为: " << fixed << setprecision(9) << res << endl;
		}

		//求阶乘的情况
		else if (Judge_factorial(init_exp))
		{
			int n = cal_factorial_number(init_exp);
			cout << "结果为: ";
			for (int i = cal_factorial(n); i > 0; i--)
				cout << LargeNumberFactorial[i];
			cout << endl;
		}

		//是加减乘除的情况
		else
		{
			cin.clear();
			cin.sync();
			while (!opt.empty())
			{
				opt.pop();
			}
			while (!val.empty())
			{
				val.pop();
			}
			del_space(init_exp); //去除空格
			string cng_exp;
			cng_exp.clear();
			change(init_exp, cng_exp); //转换为后缀表达式
			compute(cng_exp);		   //计算后缀表达式
			double stdans = val.top();
			cout << "结果为:" << (double)stdans << endl;
		}
	}
	return 0;
}

/*为每一个操作符返回一个数，数代表了优先级*/
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
}

/*判断输入的操作符是否合法*/
bool in_set(char theChar)
{
	for (int i = 0; i < 8; i++)
	{
		if (theChar == opt_set[i])
			return true;
	}
	return false;
}

/*用于去除空格并检查是否有非法字符*/
bool del_space(string &theString)
{
	string res;
	for (int i = 0; i < theString.length(); i++)
	{
		if (in_set(theString[i]) || isdigit(theString[i]))
			res += theString[i];
		else if (theString[i] == ' ')
		{
		}
		else
		{
			cout << "The enter is Wrong" << endl;
			return false;
		}
	}
	theString = res;
	return true;
}
/*把数字转换成字符串*/
string to_string(int theInt)
{
	if (theInt == 0)
		return string("0");
	bool neg = false;
	if (theInt < 0)
	{
		neg = true;
		theInt = -theInt;
	}
	string res;
	while (theInt != 0)
	{
		char c = (theInt % 10) + '0';
		res = c + res;
		theInt /= 10;
	}
	if (neg)
		res = '-' + res;
	return res;
}

/*将中缀表达式转换成后缀表达式*/
bool change(string &from, string &to)
{
	int theInt = 0;
	int state = OUT;
	char c;

	for (int i = 0; i < from.length(); i++)
	{
		c = from[i];
		if (isdigit(c))
		{
			to = to + c;
			theInt *= 10;
			theInt += c - '0';
			state = IN; //状态改为在数字内
		}
		else
		{
			if (state == IN && c == '.')
			{
				to = to + '.';
				theInt = 0;
				continue;
			}
			if (state == IN && c != '.')
			{
				to += ' ';
				theInt = 0;
			}
			if (c == '=')
				break;
			else if (c == '(')
				opt.push(c);
			else if (c == ')')
			{
				while (!opt.empty() && opt.top() != '(')
				{
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
	cout << to << endl;
	return true;
}

/*计算后缀表达式*/
bool compute(string &theExp)
{
	int theInt = 0;	 //暂存数字
	int state = OUT; //初始状态为在数字外
	char c;
	bool dot = false;
	double count = 1.0;
	for (int i = 0; i < theExp.length(); i++)
	{
		c = theExp[i];
		if (isdigit(c) || c == '.')
		{
			if (isdigit(c))
			{
				theInt *= 10;
				theInt += c - '0';
				state = IN; //状态为在数字内
				if (dot == true)
					count *= 10.0;
			}
			if (c == '.')
			{
				dot = true;
				continue;
			}
		}
		else
		{
			dot = false;
			double ans = theInt / count;
			count = 1.0;
			if (state == IN)
			{
				val.push(ans);
				theInt = 0;
			}
			double x, y;
			if (c != ' ')
			{
				x = val.top();
				val.pop();
				y = val.top();
				val.pop();

				switch (c)
				{
				case '+':
					val.push(x + y);
					break;
				case '-':
					val.push(y - x);
					break;
				case '*':
					val.push(x * y);
					break;
				case '/':
					val.push(y / x);
					break;
				default:
					cout << "未知的错误!" << endl;
				}
			}
			state = OUT;
		}
	}
	return true;
}

/*判断是否为三角函数或反三角函数并求出里面的数值*/
int Judge_rectangle(string &theStr)
{
	string str1 = "sin";
	string str2 = "cos";
	string str3 = "tan";
	string str4 = "arcsin";
	string str5 = "arccos";
	string str6 = "arctan";
	const char *show1, *show2, *show3, *show4, *show5, *show6;
	show4 = strstr(theStr.c_str(), str4.c_str());
	if (show4 != NULL)
		return 4;
	show5 = strstr(theStr.c_str(), str5.c_str());
	if (show5 != NULL)
		return 5;
	show6 = strstr(theStr.c_str(), str6.c_str());
	if (show6 != NULL)
		return 6;
	show1 = strstr(theStr.c_str(), str1.c_str());
	if (show1 != NULL)
		return 1;
	show2 = strstr(theStr.c_str(), str2.c_str());
	if (show2 != NULL)
		return 2;
	show3 = strstr(theStr.c_str(), str3.c_str());
	if (show3 != NULL)
		return 3;
	return 0;
}

/*求三角函数里包含的数字并求出对应的结果*/
double cal_rect_number(string &theStr, int theInt)
{
	int number = 0;
	for (int i = 0; theStr[i] != ')'; i++)
	{
		if (isdigit(theStr[i]))
		{
			number *= 10;
			number += theStr[i] - '0';
		}
	}
	double ans = 0.0;
	switch (theInt)
	{
	case 1:
		ans = sin(number);
		break;
	case 2:
		ans = cos(number);
		break;
	case 3:
		ans = tan(number);
		break;
	case 4:
		ans = asin(number);
		break;
	case 5:
		ans = acos(number);
		break;
	case 6:
		ans = atan(number);
		break;
	default:
		return 0;
		break;
	}
	return ans;
}

/*判断是否为求阶乘*/
bool Judge_factorial(string &theStr)
{
	int len = theStr.length();
	for (int i = 0; i < len; i++)
		if (theStr[i] == '!')
			return true;
	return false;
}

/*求阶乘所求的数字*/
int cal_factorial_number(string &theStr)
{
	int number = 0;
	for (int i = 0; theStr[i] != '!'; i++)
	{
		if (isdigit(theStr[i]))
		{
			number *= 10;
			number += theStr[i] - '0';
		}
	}
	return number;
}

/*利用数组求阶乘*/
int cal_factorial(int theInt)
{
	LargeNumberFactorial[0] = 1;
	LargeNumberFactorial[1] = 1;
	int digit = 1;
	for (int i = 1; i <= theInt; i++)
	{
		for (int j = 1; j <= digit; j++)
		{
			LargeNumberFactorial[j] *= i;
		}
		for (int j = 1; j <= digit; j++)
		{
			if (LargeNumberFactorial[j] > 10)
			{
				for (int k = 1; k <= digit; k++)
				{
					if (LargeNumberFactorial[digit] > 9)
						digit++;
					LargeNumberFactorial[k + 1] += LargeNumberFactorial[k] / 10;
					LargeNumberFactorial[k] %= 10;
				}
			}
		}
	}
	return digit;
}
