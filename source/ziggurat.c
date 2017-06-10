/* OneJoker RNG library <http://lcrocker.github.io/onejoker/randlib>
 *
 * To the extent possibile under law, Lee Daniel Crocker has waived all
 * copyright and related or neighboring rights to this work.
 * <http://creativecommons.org/publicdomain/zero/1.0/>
 *
 * This code is based on George Marsaglia's Ziggurat algorithm, with
 * various improvements from Jurgen A. Doornik and me. I believe it is
 * sufficiently dissimilar to any pre-existing implementation that there
 * should be no copyright problems, so my CC0 dedication applies.
 */

#include <stdint.h>
#include <math.h>

// On some (32-bit) machines, comparison of doubles might actually be
// faster than 64-bit integer compares, so comment out this define
#define INTEGER_COMPARE 1
#include "zigtables.h"
#include "xorshift1024star.h"

// The sizes of 128 and 256 from GM's paper are good, but I recalculated
// the constants here and built fixed tables based on them.

#define ZNR128 3.4426198558966521214
#define ZNV128 0.0099125630353364610791
#define ZER256 7.6971174701310497140
#define ZEV256 0.0039496598225815572200


static inline double ojr_to_double(uint64_t x) {
    union {
    	uint64_t r;
	double d;
    } u;
    u.r = (x & 0xFFFFFFFFFFFFFull) | 0x3FF0000000000000ull;
    return u.d - 1.0;
}


double ojr_next_exponential(void) {
    uint64_t r;
    int i;
    double x, u0, f0, f1;

    while (1) {
        r = xorshift1024_next();
        i = r & 0xFF;
        r = (r >> 8) & 0xFFFFFFFFFFFFFULL;

#ifdef INTEGER_COMPARE
        if (r < zeri[i]) {
	    u0 = ojr_to_double(r);
            return u0 * zex[i];
        }
#else
	u0 = ojr_to_double(r);

        if (u0 < zer[i]) return u0 * zex[i];
#endif
        if (0 == i) return ZER256 - log(ojr_to_double(xorshift1024_next()));

#ifdef INTEGER_COMPARE
	u0 = ojr_to_double(r);
#endif
        x = u0 * zex[i];
        f0 = exp(x - zex[i]);
        f1 = exp(x - zex[i+1]);
        if (f1 + ojr_to_double(xorshift1024_next()) * (f0 - f1) < 1.0) return x;
    }
}

double ojr_next_normal(void) {
    uint64_t r;
    int i, sign;
    double x, y, a, u0, f0, f1;

    while (1) {
        do {
            r = xorshift1024_next();
            sign = (int)r & 1;
            i = (r >> 1) & 0x7F;
            r >>= 12;
        } while (sign && 0LL == r);

#ifdef INTEGER_COMPARE
        if (r < znri[i]) {
	    a = ojr_to_double(r);
            return znx[i] * (sign ? -a : a);
        }
#else
	a = ojr_to_double(r);
        u0 = sign ? -a : a;

        if (a < znr[i]) return u0 * znx[i];
#endif
        if (0 == i) {
            do {
                x = log(ojr_to_double(xorshift1024_next())) / ZNR128;
                y = log(ojr_to_double(xorshift1024_next()));
            } while (-2.0 * y < x * x);
            return sign ? x - ZNR128 : ZNR128 - x;
        }
#ifdef INTEGER_COMPARE
	a = ojr_to_double(r);
        x = znx[i] * (sign ? -a : a);
#else
        x = u0 * znx[i];
#endif
        f0 = exp(-0.5 * (znx[i] * znx[i] - x * x));
        f1 = exp(-0.5 * (znx[i+1] * znx[i+1] - x * x));
        if (f1 + ojr_to_double(xorshift1024_next()) * (f0 - f1) < 1.0) return x;
    }
}
