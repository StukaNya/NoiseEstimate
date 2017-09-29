#ifndef DECODER_H
#define DECODER_H

#include <QObject>
#include <QVector>
#include <QDebug>
#include <QFile>
#include <QFileInfo>
#include <QtEndian>
#include <QLabel>

#include "img_pars.h"

typedef unsigned short int uInt16;
typedef unsigned char      uInt8;
class DecoderJPG2 : public QObject
{
    Q_OBJECT

public:
    static DecoderJPG2& instance() {
        static DecoderJPG2 s;
        return s;
    }

private:
    //qfile opening
    QFile inFile;
    QByteArray inStream;
    //qimage opening
    QImage inImg;
    QVector<QVector<uchar> > rgbData;
    //singleton methods
    DecoderJPG2();
    ~DecoderJPG2();
    DecoderJPG2(DecoderJPG2 const&) = delete;
    DecoderJPG2& operator= (DecoderJPG2 const&) = delete;
    //algorithm
    void openFile(const QString strDir);
    void convertToArray();
    //send message to QTextBox
    void sendMessage(const QString str);
    imgPars imSize;
signals:
    //send image and imgPars to
    void doNoiseEstimate(const QVector<QVector<uchar> >, const imgPars, int numProc);
    void doSendMessage(const QString str);

public slots:
    void runDecoding(const QString strDir, int numProc);
};

#endif // DECODER_H
