#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <QLayout>
#include <QtWidgets>
#include <QLabel>
#include <qcustomplot.h>

class MainWidget : public QWidget
{
    Q_OBJECT

public:
    static MainWidget& instance()
    {
        static MainWidget s;
        return s;
    }
    MainWidget(QWidget *parent = 0);
    ~MainWidget();
signals:
    void doSendParams(const QString dir, int numProc);
private:
    QPushButton *buttonLoad;
    QPushButton *buttonRun;
    QLabel      *labelDir;
    QLineEdit   *lineDir;
    QLabel      *labelNumProc;
    QSpinBox    *lineNumProc;
    QTextEdit   *outStream;
    QVBoxLayout *vLayout;
    QHBoxLayout *h1Layout;
    QHBoxLayout *h2Layout;
    QLabel      *imLabel;
    QCustomPlot *cPlot;
    QVector<float> xAx;
    QVector<float> yAx;
    bool flag;
    MainWidget(MainWidget* const&) = delete;
    MainWidget& operator= (MainWidget const&) = delete;

    void drawImage();
    QString strDir;
public slots:
    void receiveMessage(const QString str);
    void receiveAxes(QVector<float> xAxe, QVector<float> yAxe, int nChn, int szChn);

private slots:
    void sendParams();
    void loadImage();
};

#endif // WIDGET_H
