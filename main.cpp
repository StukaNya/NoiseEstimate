#define DEBUG_MAIN

#include "widget.h"
#include "decoder.h"
#include "noise_estimate.h"
#include <QApplication>

#ifdef DEBUG_MAIN
#include <QDebug>
#endif

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWidget widgetInst;
    DecoderJPG2& decoderInst = DecoderJPG2::instance();
    NoiseEstimate& noiseEstInst = NoiseEstimate::instance();
    //посылает адрес файла и число потоков в декодер
    QObject::connect(&widgetInst, &MainWidget::doSendParams,
                     &decoderInst, &DecoderJPG2::runDecoding);
    //запись сообщений из декодера в QTextBox
    QObject::connect(&decoderInst, &DecoderJPG2::doSendMessage,
                     &widgetInst, &MainWidget::receiveMessage);
    QObject::connect(&noiseEstInst, &NoiseEstimate::doSendMessage,
                     &widgetInst, &MainWidget::receiveMessage);
    //запуск алгоритма пономаренко
    QObject::connect(&decoderInst, &DecoderJPG2::doNoiseEstimate,
                     &noiseEstInst, &NoiseEstimate::runPonomarenkoEstimate);
    //посылает данные алгоритма и строит график
    QObject::connect(&noiseEstInst, &NoiseEstimate::doSendAxes,
                     &widgetInst, &MainWidget::receiveAxes);

    //соединить Qtext с обоими классами
    widgetInst.show();

    return a.exec();
}
