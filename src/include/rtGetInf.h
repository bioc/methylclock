/**
 * @file rtGetInf.h
 * @brief BNN neural-network clock: runtime infinity support.
 */

#ifndef RTGETINF_H
#define RTGETINF_H
#include <stddef.h>
#include "rtwtypes.h"
#include "rt_nonfinite.h"

extern real_T rtGetInf(void);
extern real32_T rtGetInfF(void);
extern real_T rtGetMinusInf(void);
extern real32_T rtGetMinusInfF(void);

#endif
