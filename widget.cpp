#define DEBUG_WIDGET

#include "widget.h"
#include <QFileDialog>

#ifdef DEBUG_WIDGET
#include <QDebug>
#endif

MainWidget::MainWidget(QWidget *parent)
    : QWidget(parent)
{
    //set title and geometry
    this->setWindowFlags(Qt::Window);
    this->setWindowTitle("JPEG2000 Converter");
    this->setGeometry(350, 100, 700, 400);

    //init buttons & labels
    buttonLoad      = new QPushButton("Load image");
    buttonRun       = new QPushButton("Run");
    labelDir        = new QLabel("Image directory");
    lineDir         = new QLineEdit();
    labelNumProc    = new QLabel("Number of threads");
    lineNumProc     = new QSpinBox();
    outStream       = new QTextEdit();
    imLabel         = new QLabel;
    cPlot           = new QCustomPlot();

    cPlot->setFixedSize(QSize(700, 400));
    lineNumProc->setRange(1,4);

    //connect buttons
    connect(buttonLoad, SIGNAL(clicked()), this, SLOT(loadImage()));
    connect(buttonRun, SIGNAL(clicked()), this, SLOT(sendParams()));

    //init layout
    vLayout = new QVBoxLayout;
    h1Layout = new QHBoxLayout;
    h2Layout = new QHBoxLayout;

    //add widgets to layout
    h1Layout->addWidget(labelDir);
    h1Layout->addWidget(lineDir);
    h1Layout->addWidget(buttonLoad);

    h2Layout->addWidget(labelNumProc);
    h2Layout->addWidget(lineNumProc);
    h2Layout->addWidget(buttonRun);

    vLayout->addWidget(cPlot);
    vLayout->addLayout(h1Layout);
    vLayout->addLayout(h2Layout);
    vLayout->addWidget(outStream);
    flag = false;
    this->setLayout(vLayout);
}

void MainWidget::loadImage() {
    strDir = QFileDialog::getOpenFileName(0, QObject::tr("Open Image"),
                                          "C:/Rust_pr/Kursach/Jpeg2/images",
                                          QObject::tr("Images (*.jpg *.jp2 *.png)"));
    this->lineDir->setText(strDir);
}

void MainWidget::sendParams() {
    strDir = this->lineDir->text();
    int numProc = this->lineNumProc->value();

#ifdef DEBUG_WIDGET
    qInfo() << "In widget class:";
    qInfo() << "Dir: " << strDir;
    qInfo() << "Num Proc: " << numProc;
#endif

    cPlot->clearPlottables();
    emit doSendParams(strDir, numProc);
}

void MainWidget::receiveMessage(const QString str)
{
    this->outStream->append(str);
}


void MainWidget::receiveAxes(QVector<float> xAxe, QVector<float> yAxe, int nChn, int szChn) {
    cPlot->legend->setVisible(true);
    cPlot->legend->setFont(QFont("Helvetica", 9));
    cPlot->legend->setRowSpacing(-3);
    xAx = xAxe;
    yAx = yAxe;
    flag = true;

    QVector<QCPScatterStyle::ScatterShape> shapes;
    shapes << QCPScatterStyle::ssDisc;
    shapes << QCPScatterStyle::ssDiamond;
    shapes << QCPScatterStyle::ssTriangle;

    QPen pen;
    for (int i = 0; i < nChn; i++) {
        cPlot->addGraph();
        QString cName;
        switch (i) {
            case 0: pen.setColor(QColor(255,0,0)); cName = "Red"; break;
            case 1: pen.setColor(QColor(0,255,0)); cName = "Green"; break;
            case 2: pen.setColor(QColor(0,0,255)); cName = "Blue"; break;
            default: pen.setColor(QColor(255,255,0)); break;
         }
        QVector<double> x(szChn), y(szChn);
        for (int k = 0; k < szChn; k++) {
            x[k] = xAx[i * szChn + k];
            y[k] = yAx[i * szChn + k];
        }
        cPlot->graph()->setData(x, y);
        cPlot->graph()->rescaleAxes(true);

        QObject::connect(cPlot->xAxis, SIGNAL(rangeChanged(QCPRange)), cPlot->xAxis2, SLOT(setRange(QCPRange)));
        QObject::connect(cPlot->yAxis, SIGNAL(rangeChanged(QCPRange)), cPlot->yAxis2, SLOT(setRange(QCPRange)));

        cPlot->graph()->setPen(pen);
        cPlot->graph()->setName(cName);
        cPlot->graph()->setLineStyle(QCPGraph::lsLine);
        cPlot->graph()->setScatterStyle(QCPScatterStyle(shapes.at(i), 10));

        cPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    }
    cPlot->rescaleAxes();

    cPlot->xAxis->setTicks(true);
    cPlot->yAxis->setTicks(true);

    cPlot->xAxis->setTickLabels(true);
    cPlot->yAxis->setTickLabels(true);
    cPlot->axisRect()->setupFullAxesBox();

    cPlot->replot();

    drawImage();
}


void MainWidget::drawImage() {
    //image
    QImage inImage;
    inImage.load(strDir);
    imLabel->setPixmap(QPixmap::fromImage(inImage));
    imLabel->setScaledContents(true);
    imLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    imLabel->show();
}


MainWidget::~MainWidget()
{

}
