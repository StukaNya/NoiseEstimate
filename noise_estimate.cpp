
#include "noise_estimate.h"
#include "../ponomarenko/framework/operations.cpp"
//#define DEBUG_NOISE_EST
//#undef _OPENMP

void NoiseEstimate::sendMessage(const QString str)
{
    emit doSendMessage(str);
}


//! Computes the delta matrix and returns the normalization factor theta
/*!
  \param *delta Delta matrix (mask for the low/high freqs in the block)
  \param w Block side
  \param T Number of low-freq coefficients, excluding DC
  \return theta Normalization factor for the matrix delta
*/
int NoiseEstimate::compute_delta(float *delta, int w, int T) {
    int theta = 0;
    for (int j = 0; j < w; j++)
        for (int i = 0; i < w; i++) {
            int value = (i + j < T && i + j != 0 ? 1 : 0);
            delta[j*w+i] = value;
            theta += value;
        }
    return theta;
}


//! Computes the set of variances computed form the high-frequency coefficients of the given blocks
/*!
  \param *VH Output set of variances
  \param **blocks_ptr List of pointers to the blocks
  \param *indices_VL Sorting indices for the blocks_ptr list (by low-freqs)
  \param w Block side
  \param T Number of low-freq coefficients, excluding DC
  \param K Number of blocks that should be used
  \return Length of the returned variances list
*/
int NoiseEstimate::compute_VH(float *VH, float **blocks_ptr, int *indices_VL, int w,
                              int T, int K) {
    int VH_count = 0;

    //#pragma omp parallel for
    for (int q = 0; q < w*w; q++) {
        int j = q / w;
        int i = q - j*w;

        if (i + j >= T) {
            float s = 0.0;
            for (int k = 0; k < K; k++) {
                float *block = blocks_ptr[indices_VL[k]];
                s += pow(block[q], 2); // q == j*w+i
            }
            VH[VH_count++] = s / K;
        }
    }
    return VH_count;
}


//! Return the optimal T parameter according to the given block side
/*!
  \param w Block side
  \return The optimal T parameter
*/
int NoiseEstimate::get_T(int w) {
    switch (w) {
    case 3: return 3;
    case 4: return 4;
    case 5: return 5;
    case 7: return 8;
    case 8: return 9;
    case 11: return 13;
    case 15: return 17;
    case 21: return 24;
    default:
        exit(-1);
    }
}

//! Reads all valid blocks (all neighbor pixels are different when the mask
//! is active) in the image
/*!
  \param *D Output list of blocks
  \param *u Input image
  \param Nx Length of a row in the image
  \param Ny Length of a column in the image
  \param w Block side
  \param num_blocks Number of blocks
  \param *mask Mask to determine if a pixel is valid or not
  \return Number of valid block copied into the output list
*/
void NoiseEstimate::read_all_valid_blocks(float *D,
                                          float *u,
                                          int Nx, int Ny,
                                          int w, unsigned num_blocks, int *mask) {
    //if (mask == NULL) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (unsigned q = 0; q < num_blocks; q++) {
            for (int j = 0; j < w; j++) {
                for (int i = 0; i < w; i++) {
                    D[q*w*w+j*w+i] = u[j*Nx+i+q];
                }
            }
        }
    //}
/*
    else {
        unsigned *valid_coords = new unsigned[num_blocks];
        int count_coords = 0;
        //
        for (int i = 0; i < Nx*Ny; i++) {
            if (mask[i] == 0)
                valid_coords[count_coords++] = i;
        }
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (unsigned q = 0; q < num_blocks; q++) {
            int addr = valid_coords[q];

            for (int j = 0; j < w; j++) {
                for (int i = 0; i < w; i++) {
                    D[q*w*w+j*w+i] = u[j*Nx+i+addr];
                }
            }
        }
        delete[] valid_coords;
    }
*/
}

//! Computes the mean of all given blocks
/*!
  \param *means Output list of means of blocks
  \param *blocks Input list of blocks to compute their means
  \param w Block side
  \param num_blocks Number of block in the input list
*/
void NoiseEstimate::compute_means(float *means, float *blocks, int w, int num_blocks) {
    float ONE_DIV_w2 = 1.0 / (w*w);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int b = 0; b < num_blocks; b++) {
        float mean = 0.0;
        for (int p = 0; p < w*w; p++) {
            mean += blocks[b*w*w+p];
        }
        mean *= ONE_DIV_w2;
        means[b] = mean;
    }
}

int NoiseEstimate::get_max(int *data, int N) {
    int max = data[0];
    for (int i = 1; i < N; i++)
        if (data[i] > max)
            max = data[i];
    return max;
}

void NoiseEstimate::copy(float *dest, float *orig, int N) {
    for (int i = 0; i < N; i++)
        dest[i] = orig[i];
}

//! Determines if the given string corresponds to the custom percentile code
/*!
  \param *num_str Input string
  \return true if the input string corresponds to the custom percentile code or false if not.
*/
bool NoiseEstimate::is_custom_percentile(const char *num_str) {
    char buffer[1024];
    float value = atof(num_str);
    sprintf(buffer, "%.4f", value);
    return strcmp(buffer, "0.0000") != 0;
}

//! In-place Normalization of the FFTW output in order to get a orthonormal 2D DCT-II
/*!
  \param *blocks Input/output list of transformed blocks
  \param w Block side
  \param num_blocks Number of blocks in the list
*/
void NoiseEstimate::normalize_FFTW_bak(float *blocks, int w, int num_blocks) {
    const float ONE_DIV_2w = 1.0 / (2.0*w);
    const float ONE_DIV_SQRT_2 = 1.0 / sqrtf(2);

    // Divide all coefficients by 2*w
    //#pragma omp parallel for shared(blocks)
    for (int i = 0; i < num_blocks*w*w; i++)
        blocks[i] *= ONE_DIV_2w;

#ifdef _OPENMP
#pragma omp parallel for shared(blocks) schedule(static)
#endif
    for (int b = 0; b < num_blocks; b++) {
        // {(i, j)} with i == 0 or j == 0
        for (int j = 0; j < w; j++) {
            int i = 0;
            blocks[b*w*w+j*w+i] *= ONE_DIV_SQRT_2;
        }
        for (int i = 0; i < w; i++) {
            int j = 0;
            blocks[b*w*w+j*w+i] *= ONE_DIV_SQRT_2;
        }
    }
}

/**
 * @brief Build a mask for valide pixel. If mask(i, j) = true, the pixels will not be used.
 *
 * @param i_im : noisy image;
 * @param o_mask : will contain the mask for all pixel in the image size;
 * @param p_imSize : size of the image;
 * @param p_sizePatch : size of a patch.
 *
 * @return number of valid blocks.
 *
 **/
unsigned NoiseEstimate::build_mask(QVector<QVector<float> > &i_im, int *o_mask,
                                    unsigned Nx, unsigned Ny, unsigned w,
                                    unsigned num_channels) {
    unsigned count  = 0;

    for (unsigned ij = 0; ij < Nx*Ny; ij++) {
        const unsigned j = ij / Nx;
        const unsigned i = ij - j * Nx;

        //! Look if the pixel is not to close to the boundaries of the image
        if (i < Nx - w + 1 && j < Ny - w + 1) {
            for (unsigned c = 0; c < num_channels; c++) {
                float *u = &i_im[c][0];

                //! Look if the square 2x2 of pixels is constant
                int invalid_pixel = (c == 0 ? 1 : o_mask[ij]);

                // Try to validate pixel
                if (fabs(u[ij] - u[ij + 1]) > 0.001f) {
                    invalid_pixel = 0;
                }
                else
                    if (fabs(u[ij + 1] - u[ij + Nx]) > 0.001f) {
                        invalid_pixel = 0;
                    }
                    else
                        if (fabs(u[ij + Nx] - u[ij + Nx + 1]) > 0.001f) {
                            invalid_pixel = 0;
                        }
                o_mask[ij] = invalid_pixel;
            }
        }
        else {
            o_mask[ij] = 1; // Not valid
        }

        if (o_mask[ij] == 0)
            count++;
    }

    return count;
}


/**
 * @brief Ponomarenko-diff noise estimation algorithm.
 *
 * @param i_im : input noisy image (RGB format);
 * @param p_imSize : size of the image;
 * @param io_nlbParams : see nlbParams. The covariance matrices will be updated.
 * @param p_msdParams : see MsdParams.
 *
 * @return EXIT_FAILURE in case of problems, EXIT_SUCCESS otherwise.
 **/
int NoiseEstimate::runPonomarenkoEstimate(
        QVector<QVector<uchar> > i_im,
        const imgPars& p_imSize, int numProc) {

    //! Parameters initialization
    algParams params;
    params.sizePatch        = 4;
    params.percentile       = 0.005f;
    params.correctionFactor = 0.5f * 7.460923461f; // Work only for sP = 4. Needs to be tabulated for other values.
    params.numBins          = std::max(p_imSize.wh / 42000.f, 1.f);
    const int meanMethod    = 2;
    int i = 0;

    //! For convenience
    const unsigned int sP  = params.sizePatch;
    const unsigned int sP2 = sP * sP;

    // Read parameters
    int w = params.sizePatch;
    int T = (w < 6) ? w : ((w < 10) ? w - 1 : w - 2);
    float p = 0.005f;


    // Get image properties
    int Nx = p_imSize.width;
    int Ny = p_imSize.height;
    int num_channels = p_imSize.nChn;
    int total_blocks = (Nx-w+1) * (Ny-w+1); // Number of overlapping blocks

    // Parallelization config
#ifdef _OPENMP
    double start_time, stop_time;
    omp_set_num_threads(numProc);
#endif

    QVector<QVector<float> > im;
    im.resize(p_imSize.nChn);

    for (unsigned int c = 0; c < p_imSize.nChn; c++) {
        im[c].resize(p_imSize.wh);
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (unsigned int k = 0; k < p_imSize.wh; k++) {
            im[c][k] = (float)i_im[c][k];
        }
    }

    // Arrays to store the final means and noise estimations
    int *mask_all   = nullptr;
    float *delta    = nullptr;
    float *vmeans   = nullptr;
    float *vstds    = nullptr;
    //float *MatStds  = nullptr;
    float *means    = nullptr;
    float *blocks   = nullptr;

    try {
        //mask_all    = new int[Nx*Ny];
        delta       = new float[w*w];
        vmeans      = new float[num_channels * params.numBins];
        vstds       = new float[num_channels * params.numBins];
        //MatStds     = new float[num_channels * params.numBins * w*w];
    }
    catch (...) {
#ifdef DEBUG_NOISE_EST
        qInfo() << "Bad alloc! (mask_all, delta, ...)";
#endif
        //delete[] mask_all; mask_all = nullptr;
        delete[] delta  ; delta   = nullptr;
        delete[] vmeans ; vmeans  = nullptr;
        delete[] vstds  ; vstds   = nullptr;
        //delete[] MatStds; MatStds = nullptr;
    }

    int num_blocks = Nx*Ny;//build_mask(im, mask_all, Nx, Ny, w, num_channels);
    int theta = compute_delta(delta, w, T);

    try {
        means       = new float[num_blocks];
        blocks      = new float[num_blocks * w*w];
    }
    catch (...) {
#ifdef DEBUG_NOISE_EST
        qInfo() << "Bad alloc! (means, blocks)";
#endif
        delete[] means ; means  = nullptr;
        delete[] blocks; blocks = nullptr;
    }

    // Init FFTW threads
    fftwf_init_threads();

    int nbTable[2] = {w,w};
    int nembed[2] = {w,w};

#ifdef _OPENMP
    fftwf_plan_with_nthreads(numProc);
#endif

    fftwf_r2r_kind kindTable[2] = {FFTW_REDFT10, FFTW_REDFT10};

    fftwf_plan fft_plan = fftwf_plan_many_r2r(2, nbTable, num_blocks, blocks,
                                              nembed, 1, w*w, blocks, nembed,
                                              1, w*w, kindTable, FFTW_ESTIMATE);

#ifdef DEBUG_NOISE_EST
    qInfo() << "Run Ponomarenko estimate! Channels: " << num_channels;
#endif
#ifdef _OPENMP
    start_time = omp_get_wtime();
#endif
    // Process each channel
    for (int ch = 0; ch < num_channels; ch++) {
        float *u = (im[ch]).data();
        // Create a list of pointers of the groups
        float **blocks_ptr;
        try {
            blocks_ptr = new float*[num_blocks];
        }
        catch (...) {
#ifdef DEBUG_NOISE_EST
            qInfo() << "Bad alloc! In channel: " << ch;
#endif
            delete[] blocks_ptr; blocks_ptr = nullptr;
            break;
        }

        for (int i = 0; i < num_blocks; i++)
            blocks_ptr[i] = &blocks[i*w*w];

        read_all_valid_blocks(blocks, u, Nx, Ny, w, num_blocks, mask_all);

        // Compute means
        //compute_means(means, blocks, w, num_blocks);

        // Compute 2D-DCT of all the blocks
        //
        // Transform blocks with FFTW
        fftwf_execute_r2r(fft_plan, blocks, blocks);

        // Normalize FFTW output
        normalize_FFTW_bak(blocks, w, num_blocks);

        HistMean<unsigned char, unsigned int> histInst((unsigned int)Nx, (unsigned int)Ny);

        histInst.agcHistMean((i_im[ch]).data());
        histInst.sumHistMean();
        histInst.sortHistMean();

#ifdef DEBUG_NOISE_EST
        qInfo() << "Channel:" << ch << "; Num bins:" << histInst.histP.numBins;
#endif

            // Process each bin
#ifdef _OPENMP
      #pragma omp parallel for shared(vmeans, vstds, histInst) schedule(dynamic)
#endif
        for (int bin = 0; bin < histInst.histP.numBins; bin++) {
            int elems_bin = histInst.histP.samplesPerBin[bin];
            int K = elems_bin * p;
            float bin_mean, tilde_sigma;

            float *VL           = nullptr;
            float *VH           = nullptr;
            int *indices_VL     = nullptr;

            try {
                VL          = new float [elems_bin];
                indices_VL  = new int   [elems_bin];
                VH          = new float [w*w];
            }
            catch (...) {
#ifdef DEBUG_NOISE_EST
                qInfo() << "Bad alloc! In bin: " << bin;
#endif
                delete[] VL;         VL         = nullptr;
                delete[] indices_VL; indices_VL = nullptr;
                delete[] VH;         VH         = nullptr;
                //break;
            }
            //#pragma omp critical
            //{
            histInst.computeVLhist(VL, bin, w, delta, blocks_ptr, theta);

            // Compute VH
            argsort(VL, indices_VL, elems_bin);
            int VH_count = histInst.computeVHhist(VH, blocks_ptr, bin, indices_VL, w, T, K);

            //     float bin_mean = get_bin_mean(meanMethod, K, indices_VL, bin, &histo);
            bin_mean = (float)(histInst.arrMean[histInst.arrSort[bin*elems_bin+elems_bin/2]]) / 256.f;
            tilde_sigma = sqrt(median(VH, VH_count));

#ifdef DEBUG_NOISE_EST
            qInfo() << "bin:" << bin << "; mean:" << bin_mean << "; sigma:" << tilde_sigma;
#endif

            // Store results
            vmeans[ch*params.numBins+bin] = bin_mean;
            vstds[ch*params.numBins+bin] = tilde_sigma;
            //}
           // for (int ij = 0; ij < w*w; ij++) {
           //     MatStds[ch*w*w*params.numBins + bin*w*w +ij] = tilde_sigma;
           // }

            delete[] VL;         VL         = nullptr;
            delete[] VH;         indices_VL = nullptr;
            delete[] indices_VL; VH         = nullptr;
        }

        delete[] blocks_ptr; blocks_ptr = nullptr;
    }


#ifdef _OPENMP
    stop_time = omp_get_wtime();
#ifdef DEBUG_NOISE_EST
    qInfo() << "Algorithm done! Total time: " << stop_time - start_time;
#endif
    sendMessage("Total time of algorithm: " + QString::number(stop_time - start_time));
#endif



    QVector<float> xAxe(params.numBins * num_channels),
                    yAxe(params.numBins * num_channels);
    qCopy(vmeans, vmeans+params.numBins*num_channels, xAxe.begin());
    qCopy(vstds, vstds + params.numBins*num_channels, yAxe.begin());
    emit doSendAxes(xAxe, yAxe, num_channels, params.numBins);

    // FFTW Cleanup
    fftwf_destroy_plan(fft_plan);
    fftwf_cleanup_threads();
    fftwf_cleanup();

    // Clean up memory
    delete[] mask_all; mask_all = nullptr;
    delete[] delta   ; delta    = nullptr;
    delete[] vmeans  ; vmeans   = nullptr;
    delete[] vstds   ; vstds    = nullptr;
    //delete[] MatStds ; MatStds  = nullptr;
    //delete[] means   ; means    = nullptr;
    delete[] blocks  ; blocks   = nullptr;

    return 0;
}
