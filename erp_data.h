#ifndef ERP_DATA_H
#define ERP_DATA_H

#include <QVector>

class erp
{
public:
    erp();
    double MJD;
    double Xpole;
    double Ypole;
    double UT1_UTC;
    double dpsi;
    double deps;
};

class erp_file
{
public:
    erp_file();
    QVector<erp> record;
};

#endif // ERP_DATA_H
