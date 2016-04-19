#ifndef READFILEPATH_H
#define READFILEPATH_H

#include <QString>
#include <QVector>

class readFilePath
{
public:
    readFilePath();
    ~readFilePath();

    static void readFile(QVector<QString> &path , QString filePathName ,QString fileType);
};

#endif // READFILEPATH_H
