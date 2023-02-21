#include <stdio.h>
#include <iostream>
#include <float.h>
using namespace std;

















int main()
{
    cout.precision(100);
    float a = FLT_MAX;
    cout << ((a+100.0f)>FLT_MAX) << endl;
    cout << a*0.75f + a << endl;
    cout << ((a*0.75f+a)>FLT_MAX) << endl;
    return 0;
}