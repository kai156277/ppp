#include<QString>
#include<QVector>
#include<QFile>
#include<QTextStream>
#include<qdebug.h>

#include "output.h"

output::output()
{

}

void output::writePPP(const result_file &file)
{
    QFile result_file( "C:/PPP/readFile/ppp.r" );
    if(!result_file.open(QIODevice::WriteOnly))
    {
        qDebug() << "Can`t open" << "C:/PPP/readFile/ppp.r" << endl;
        exit(EXIT_FAILURE);
    }
    QTextStream write( &result_file );
    for(int i = 0; i<file.file.size(); i++)
    {
        write << file.file[i].dx << "\t"
              << file.file[i].dy << "\t"
              << file.file[i].dz << endl;
    }
}

