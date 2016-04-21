#include<iostream>
#include<QVector>

#include"file_read.h"
#include"o_date.h"
#include"sp3_date.h"
#include"clock_date.h"
#include"math_function.h"
#include"readfilepath.h"
#include"ppp_calculate.h"
#include"ppp_date.h"
using namespace std;

int main()
{
    o_file_date ofile;
    sp3_file sp3;
    file_read read;
    clock_file clock;
    ppp_calculate pppCalculate;
    ppp_file ppp;
    double d[8] ={-1,0,1,2,-7,-4,5,26};
    /*
     * QVector<QString> filePath;
     * QString pathName ="C://PPP//readFile//sp3";
     * QString fileType = "*.sp3";
     * readFilePath::readFile(filePath,pathName,fileType);
     */
    read.ppp_o_read("C://PPP//readFile//oFile//cas11680.15o",ofile);
    read.ppp_sp3_read("C:/PPP/readFile/sp3/gbm18493.sp3",sp3);
    read.ppp_clock_read("C:/PPP/readFile/clk/gbm18493.clk",clock);
    pppCalculate.ppp_clock(ofile,sp3,clock,ppp);
    //pppCalculate.ppp_coordinate(ofile,sp3,ppp);
    //double i = math_function::lagrange(0.5,d,8);
    int j = 0;
    cout << "END!";
}
