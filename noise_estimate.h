#ifndef _PONOMARENKO_H_
#define _PONOMARENKO_H_

#include <QObject>
#include <QVector>
#include <QDebug>
#include <QtAlgorithms>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "img_pars.h"
#include "curve_filter.h"
#include "hist_mean.h"
#include "decoder.h"

class NoiseEstimate : public QObject
{
    Q_OBJECT

public:
    static NoiseEstimate& instance() {
        static NoiseEstimate s;
        return s;
    }

    typedef struct algParams_t{
      float sizePatch;
      float percentile;
      uint numBins;
      uint numBlocks;
      uint currentChannel;
      float correctionFactor;
    } algParams;

private:
    int compute_delta(float *delta, int w, int T);
    int compute_VH(float *VH, float **blocks_ptr, int *indices_VL, int w, int T, int K);
    int get_T(int w);
    void read_all_valid_blocks(float *D, float *u, int Nx, int Ny, int w, unsigned num_blocks, int *mask);
    void compute_means(float *means, float *blocks, int w, int num_blocks);
    int get_max(int *data, int N);
    void copy(float *dest, float *orig, int N);
    bool is_custom_percentile(const char *num_str);
    void normalize_FFTW_bak(float *blocks, int w, int num_blocks);
    unsigned build_mask(QVector<QVector<float> > &i_im, int *o_mask,
                        unsigned Nx, unsigned Ny, unsigned w,
                        unsigned num_channels);
    void sendMessage(const QString str);
signals:
    void doSendMessage(const QString str);
    void doSendAxes(QVector<float> xAxe, QVector<float> yAxe, int nChn, int szChn);

public slots:
    int runPonomarenkoEstimate(
        QVector<QVector<uchar> > i_im,
        const imgPars& p_imSize, int numProc);
};
#endif


