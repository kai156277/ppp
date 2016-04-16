#include<iostream>
#include<QVector>

#include"o_date.h"
#include"file_read.h"
#include"sp3_date.h"
#include"clock_date.h"
#include"sate_c_c.h"
#include"math_function.h"
using namespace std;

int main()
{
    o_file_date a;
    sp3_file b;
    file_read test;
    clock_date c;
    double d[8] = {-1,0,1,2,-7,-4,5,26};
    //test.ppp_o_read("C://PPP//readFile//oFile//cas11680.15o",a);
    //test.ppp_sp3_read("C:/PPP/readFile/sp3/gbm18492.sp3",b);
    //test.ppp_clock_read("C:/PPP/readFile/clk/gbm18492.clk",c);
    double i = math_function::lagrange(0.5,d,4);
    int j = 0;
}
