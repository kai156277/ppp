#ifndef PPP_CALCULATE_H
#define PPP_CALCULATE_H

#include<QVector>
#include<QString>

#include"o_date.h"
#include"sp3_date.h"
#include"clock_date.h"
#include"ppp_date.h"

class ppp_calculate
{
public:
    ppp_calculate();
    void ppp_coordinate(const o_file_date &ofile,const sp3_file &sp3file,ppp_file &ppp);
    void ppp_clock(const o_file_date &ofile,const sp3_file &sp3file,const clock_file &clockfile,ppp_file &ppp);
private:
    static double c;
};

#endif // PPP_CALCULATE_H
