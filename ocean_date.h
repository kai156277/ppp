#ifndef OCEAN_DATE_H
#define OCEAN_DATE_H

#include<QStringList>
#include<QVector>


class ocean
{
public:
    ocean();
    QString stationName;
    QString oceanModel;
    int year;
    int month;
    int day;
    double station_B;
    double station_L;
    double station_H;
    double data[6][11];
};

class ocean_file
{
public:
    ocean_file();
    int col;
    QStringList signer;
    QVector<ocean> record;
};


#endif // OCEAN_DATE_H
