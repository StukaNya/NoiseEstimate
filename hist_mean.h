/*
Гистограмма средних значений (для блоков 4х4 пикселя)
Содержит:
- массив гистограммы, среднее округленное (до целого) значение равно индексу
- массив кумулятивной суммы (+ функция заполнения массива)
- массив отсортированных значений индексов (+функция сортировки)
- функция заполнения гистограммы картинкой (интенсивности: 0x0-0x3fff)
*/
#ifndef _HIST_MEAN_H_
#define _HIST_MEAN_H_

#define SIZE_QUEUE		0x3fff	//2^14
#define WIDTH			640
#define HEIGH			512
#define PATCH_SIZE		4
#define IMG_SIZE		(WIDTH-PATCH_SIZE)*(HEIGH-PATCH_SIZE)

template <typename dataT, typename argT>
class HistMean {

typedef struct histPars_t
{
    //image size
    unsigned int imgWidth;
    unsigned int imgHeight;
    unsigned int imgSz;
    //hist size
    unsigned int histSz;
    //number of bins
    unsigned int numBins;
    //кол-во элементов внутри блока
    unsigned int elemsBin;
	//количество элементов в каждом блоке
    unsigned int *samplesPerBin;
} histPars;

public:
    HistMean(unsigned int img_width, unsigned int img_height);
    ~HistMean();
    //заполнение гистограммы
    void agcHistMean(dataT *img);
    //вычисление кумулятивной суммы
    void sumHistMean();
    //сортировка, возвращает массив отсортированных индексов
    void sortHistMean();
    //VL и VH для новой гистограммы
    void computeVLhist(float *VL, int bin, int w, float *delta, float **blocks_ptr, int theta);
    int computeVHhist(float *VH, float **blocks_ptr, int bin, int *indices_VL, int w, int T, int K);

    //instance of parameters struct
    histPars histP;
    //array of means
    dataT *arrMean;
    //array of sorted indexes
    argT *arrSort;
private:
    //pointer to the indices of each block
    argT **binPtr;
    //array of histogram
    argT *arrIdx;
    //array of cumulative sum
    argT *arrSum;
};


template <typename dataT, typename argT>
HistMean<dataT, argT>::HistMean(unsigned int img_width, unsigned int img_height) {
    int divNumBin;
    int i;
    int tempIdx = 0;

    //init params
    histP.imgWidth = img_width;
    histP.imgHeight = img_height;
    histP.imgSz = img_width * img_height;
    histP.histSz = 1 << (sizeof(dataT) * 8);
    histP.numBins = (histP.imgSz / 42000.f > 1.f) ? (int)(histP.imgSz / 42000.f) : 1;
    histP.samplesPerBin = new unsigned int[histP.numBins];

    //allocate arrays
    try {
        binPtr  = new argT*[histP.numBins];
        arrMean = new dataT[histP.imgSz];
        arrIdx  = new argT[histP.histSz];
        arrSum  = new argT[histP.histSz];
        arrSort = new argT[histP.imgSz];
    }
    catch (...) {
        delete[] binPtr ; binPtr  = nullptr;
        delete[] arrMean; arrMean = nullptr;
        delete[] arrIdx ; arrIdx  = nullptr;
        delete[] arrSum ; arrSum  = nullptr;
        delete[] arrSort; arrSort = nullptr;
    }

    divNumBin = (int) (histP.imgSz / histP.numBins);
    for (i = 0; i < histP.numBins; i++)
    {
        histP.samplesPerBin[i] = (i < histP.numBins - 1) ? divNumBin : histP.imgSz - (histP.numBins - 1) * divNumBin;
        //каждый указатель выделяет блок отсортированных индексов длины samplesPerBin[i]
        binPtr[i] = &(arrSort[tempIdx]);
        tempIdx += histP.samplesPerBin[i];
    }
}


//заполнение массива кумулятивных сумм
template <typename dataT, typename argT>
void HistMean<dataT, argT>::sumHistMean() {
    argT i, sum = 0;
    arrSum[0] = 0;
    for (i = 0; i < histP.histSz - 1; i++)
    {
        sum += arrIdx[i];
        arrSum[i+1] = sum;

    }
}


//сортировка при помощи гистограммы (за один проход по изображению)
template <typename dataT, typename argT>
void HistMean<dataT, argT>::sortHistMean() {
    // i и idx - индексы пикселя до и после сортировки
    argT i, idx;
    dataT mean;
    for (i = 0; i < histP.imgSz; i++)
    {
        arrSort[i] = 0;
    }
    for (i = 0; i < histP.imgSz; i++)
    {
        mean = arrMean[i];
        idx = arrSum[mean]++;
        arrSort[idx] = i;
    }
}


//заполнение гистограммы (среднее значение блока округлено до целого)
template <typename dataT, typename argT>
void HistMean<dataT, argT>::agcHistMean(dataT *img) {

    int n, m, i, j;
    unsigned int width = histP.imgWidth;
    unsigned int height = histP.imgHeight;
    const int shiftW = 2;
    const int shiftW_2 = 2 * shiftW;
    const int w = 1 << shiftW;
    unsigned int tempSum = 0;
    dataT mean;

    for (n = 0; n < histP.histSz; n++)
    {
        arrIdx[n] = 0;
        arrSum[n] = 0;
    }

    for (n = 0; n < width - w; n++)
    {
        for (m = 0; m < height - w; m++)
        {
            tempSum = 0;
            for (i = 0; i < w; i++)
            {
                for (j = 0; j < w; j++)
                {
                    tempSum += img[(m + j) * width + n + i];
                }
            }
            mean = (dataT)(tempSum >> shiftW_2);
            arrIdx[mean]++;
            arrMean[m*width + n] = mean;
        }
    }
}


template <typename dataT, typename argT>
HistMean<dataT, argT>::~HistMean() {
    delete[] binPtr;  binPtr = nullptr;
    delete[] arrMean; arrMean= nullptr;
    delete[] arrIdx;  arrIdx = nullptr;
    delete[] arrSum;  arrSum = nullptr;
    delete[] arrSort; arrSort= nullptr;
}


//compute VL coeffs with new histogram
template <typename dataT, typename argT>
void HistMean<dataT, argT>::computeVLhist(float *VL, int bin, int w, float *delta, float **blocks_ptr, int theta) {
    int i, j, m;
    for (m = 0; m < histP.samplesPerBin[bin]; m++)
    {
        int idx = binPtr[bin][m];
        float *block = blocks_ptr[idx];
        VL[m] = 0;
        for (i = 0; i < w; i++)
        {
            for (j = 0; j < w; j++)
                if (delta[j*w+i] != 0)
                    VL[m] += block[j*w+i]*block[j*w+i];
        }
        VL[m] /= theta;
    }
}


template <typename dataT, typename argT>
int HistMean<dataT, argT>::computeVHhist(float *VH, float **blocks_ptr, int bin, int *indices_VL, int w,
                              int T, int K) {
    int VH_count = 0;

    //#pragma omp parallel for
    for (int q = 0; q < w*w; q++) {
        int j = q / w;
        int i = q - j*w;

        if (i + j >= T) {
            float s = 0.0;
            for (int k = 0; k < K; k++) {
                float *block = blocks_ptr[binPtr[bin][indices_VL[k]]];
                s += pow(block[q], 2); // q == j*w+i
            }
            VH[VH_count++] = s / K;
        }
    }
    return VH_count;
}
#endif
