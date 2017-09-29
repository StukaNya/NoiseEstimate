/*
Гистограмма средних значений (для блоков 4х4 пикселя)
Содержит:
- массив гистограммы, среднее округленное (до целого) значение равно индексу
- массив кумулятивной суммы (+ функция заполнения массива)
- массив отсортированных значений индексов (+функция сортировки)
- функция заполнения гистограммы картинкой (интенсивности: 0x0-0x3fff)
*/

#include "hist_mean.h"

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
    binPtr = new argT*[histP.numBins];
    arrMean = new dataT[histP.imgSz];
    arrIdx  = new argT[histP.histSz];
    arrSum  = new argT[histP.histSz];
    arrSort = new argT[histP.imgSz];

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
    delete[] binPtr;
    delete[] arrMean;
    delete[] arrIdx;
    delete[] arrSum;
    delete[] arrSort;
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


