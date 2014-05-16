#ifndef __COLOR_VECTOR__
#define __COLOR_VECTOR__

#include "amplitudes.h"

typedef struct { double n,d;  } rat;

extern void generateColorVectors (vamplExt * ans, int * nvect, rat ** c_vectors);

#endif
