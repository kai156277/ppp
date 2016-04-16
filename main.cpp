#include<iostream>

#include"o_date.h"
#include"file_read.h"
#include"sp3_date.h"
#include"clock_date.h"
using namespace std;

int main()
{
    o_file_date a;
    sp3_file b;
    file_read test;
    clock_date c;
    //test.ppp_o_read("C://PPP//readFile//oFile//cas11680.15o",a);
    //test.ppp_sp3_read("C:/PPP/readFile/sp3/gbm18492.sp3",b);
    test.ppp_clock_read("C:/PPP/readFile/clk/gbm18492.clk",c);
}
