#ifndef OUTPUT_H
#define OUTPUT_H

#include"ppp_data.h"

class output
{
public:
    output();
    static void writePPP(const result_file &file);
};

#endif // OUTPUT_H
