#ifndef READFILEPATH_H
#define READFILEPATH_H

#include <QString>
#include <QVector>

class readFilePath
{
public:
    readFilePath();
    ~readFilePath();

    void readObsFilePath(QVector<QString> &path , QString filePathName);
};

#endif // READFILEPATH_H
