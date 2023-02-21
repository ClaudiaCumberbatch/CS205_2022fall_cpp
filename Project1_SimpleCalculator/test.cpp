#include <iostream>
#include <string>
#include <cmath>
#include <vector>
using namespace std;

int main(int argc, char * argv[]) 
{
    long long a1 = 12345678900ll;
    long long result =(long long) (a1 * a1);
    cout.precision(50);
    cout << a1 << endl;
    cout << a1 * a1 << endl;
}

// #include <iostream>
// #include <sstream>
// #include <iomanip>

// int main(int argc, char *argv[])
// {
//     double d = 3.1415926535897932384;
//     std::string str = std::to_string(d);
//     std::cout << str << std::endl; // 3.141593

//     std::stringstream ss;
//     ss << std::setprecision(15) << d;
//     str = ss.str(); 
//     std::cout << str << std::endl; //3.14159265358979

//     return 0;
// }