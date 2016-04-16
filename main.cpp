#include<iostream>

#include"o_date.h"
#include"file_read.h"
#include"sp3_date.h"
using namespace std;

int main()
{
    o_file_date a;
    sp3_file b;
    file_read test;
    test.ppp_o_read("C://PPP//readFile//oFile//cas11680.15o",a);
    test.ppp_sp3_read("C:/PPP/readFile/sp3/gbm18492.sp3",b);
}
