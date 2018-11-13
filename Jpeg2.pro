#-------------------------------------------------
#
# Project created by QtCreator 2017-04-07T12:40:39
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Jpeg2
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += main.cpp\
        widget.cpp \
    decoder.cpp \
    noise_estimate.cpp \
    qcustomplot.cpp

HEADERS  += widget.h \
    decoder.h \
    noise_estimate.h \
    hist_mean.h \
    img_pars.h \
    qcustomplot.h

RESOURCES += \

QMAKE_LIBS += -lgomp
QMAKE_CXXFLAGS += -fopenmp -std=c++11
QMAKE_LFLAGS += -fopenmp

INCLUDEPATH += C:/Qt/Qt5.8.0/library/fftw
LIBS += -LC:/Qt/Qt5.8.0/library/fftw -lfftw3 \
        -LC:/Qt/Qt5.8.0/library/fftw -lfftw3f \
        -LC:/Qt/Qt5.8.0/library/fftw -lfftw3l

