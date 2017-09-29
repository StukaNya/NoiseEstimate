#ifndef IMG_PARS_H
#define IMG_PARS_H

#include <QMetaType>

struct imgPars {
    uint wh;
    uint width;
    uint height;
    uint nChn;
    uint bitsPChn;
};

Q_DECLARE_METATYPE(imgPars)

#endif // IMG_PARS_H
