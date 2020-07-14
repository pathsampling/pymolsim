#include <sstream>
#include <string>

using namespace std;
string num2str(double num){
    stringstream ss;
    ss<<num;
    return ss.str();
}
