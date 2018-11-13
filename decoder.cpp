#define DEBUG_DECODER

#include "decoder.h"


#ifdef DEBUG_DECODER
#include <QDebug>
#endif

DecoderJPG2::DecoderJPG2()
{

}

void DecoderJPG2::sendMessage(const QString str)
{
    emit doSendMessage(str);
}

void DecoderJPG2::runDecoding(const QString strDir, int numProc)
{
    const QString message = "Dir: " + strDir + "; Num proc: " + QString::number(numProc);
    sendMessage(message);

#ifdef DEBUG_DECODER
    qInfo() << "In decoder class:";
    qInfo() << "Dir: " << strDir;
    qInfo() << "Num Proc: " << numProc;
#endif
    if (strDir != currentDir) {
        //openFile(strDir);
        inImg.load(strDir);
        if (inImg.format() != QImage::Format_RGB32) {
            inImg = inImg.convertToFormat(QImage::Format_RGB32);
        }

        imSize.width   = inImg.width();
        imSize.height  = inImg.height();
        imSize.wh      = imSize.width * imSize.height;
        imSize.nChn    = (inImg.format() == QImage::Format_RGB32) ? 3 : 1;

    #ifdef DEBUG_DECODER
        qInfo() << "inImage size: " << inImg.size();
        qInfo() << "inImage format: "<< inImg.format();
    #endif

        convertToArray();
        currentDir = strDir;
    }
    emit doNoiseEstimate(rgbData, imSize, numProc);
}


void DecoderJPG2::convertToArray() {
    //uInt8 *blockPtr = (uInt8*)(inStream.data());
    //jp2h - header
    //ihdr - image header, 0x012C = 300 - width, 0x0190 = 400 - height
   int i;
   uint *iter = (uint *)(inImg.bits());

   rgbData.resize(imSize.nChn);
   for (i = 0; i < imSize.nChn; i++)
       rgbData[i].resize(imSize.wh);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
   for (i = 0; i < imSize.wh; i++)
   {
       rgbData[0][i] = (uchar)(*iter);
       rgbData[1][i] = (uchar)((*iter) >> 8);
       rgbData[2][i] = (uchar)((*iter) >> 16);
       iter++;
   }
}


DecoderJPG2::~DecoderJPG2()
{
    if(inFile.isOpen()) {
        inFile.close();
    }
}


void DecoderJPG2::openFile(const QString strDir)
{
    inFile.setFileName(strDir);
    if (!inFile.exists()) {
        sendMessage("File does not exist!");
        return;
    }

    QString errMsg;
    if (!inFile.open(QIODevice::ReadOnly)) {
        errMsg = inFile.errorString();
        sendMessage(errMsg);
        return;
    }

    inStream = inFile.readAll();

#ifdef DEBUG_DECODER
    qInfo() << "Img size in bytes: " << inStream.size();
#endif
}
