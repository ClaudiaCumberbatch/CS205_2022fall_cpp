#include <cctype>
#include <iostream>
#include <cstring>
using namespace std;

int main()
{
  int * p = new int[5]();
  p[0] = 1;
  cout << p[0] << endl;
  cout << p << endl;
  delete []p;
  p = NULL;

  p = new int[5]();
  p[0] = 2;
  cout << p[0] << endl;
  cout << p << endl;

  return 0;

  return 0;
}