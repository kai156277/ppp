#include<iostream>
#include<QVector>

#include"file_read.h"
#include"o_data.h"
#include"sp3_data.h"
#include"clock_data.h"
#include"math_function.h"
#include"readfilepath.h"
#include"ppp_calculate.h"
#include"ppp_data.h"
#include"snx_data.h"
#include"antmod_data.h"
using namespace std;

int main()
{
    o_file ofile;
    sp3_file sp3;
    file_read read;
    clock_file clock;
    ppp_calculate pppCalculate;
    ppp_file ppp;
    snx_data snx;
    antmod_file ant;
    /*
     * QVector<QString> filePath;
     * QString pathName ="C://PPP//readFile//sp3";
     * QString fileType = "*.sp3";
     * readFilePath::readFile(filePath,pathName,fileType);
     */
    read.ppp_o_read("C://PPP//readFile//oFile//cut01680.15o",ofile);
    read.ppp_sp3_read("C:/PPP/readFile/sp3/gbm18493.sp3",sp3);
    read.ppp_clock_read("C:/PPP/readFile/clk/gbm18493.clk",clock);
    read.ppp_ant_read("C:/PPP/readFile/ant/antmod.dat",ant);

    //read.ppp_snx_read("C:/PPP/readFile/snx/igs14P1783.snx",snx,ofile.heard.marker_name);
    pppCalculate.ppp_pretreatment(ofile,ant);


    pppCalculate.ppp_spp(ofile,sp3,clock,ant,ppp);
    /**/

    cout << "END!";
}
