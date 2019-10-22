/*  Written in 2014 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>
#include <string.h>
#include "splitmix64.h"

/* This is a fast, top-quality generator. If 1024 bits of state are too
   much, try a xoroshiro128+ generator.

   Note that the three lowest bits of this generator are LSFRs, and thus
   they are slightly less random than the other bits. We suggest to use a
   sign test to extract a random Boolean value.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static uint64_t xorshift1024_s[16];
static int xorshift1024_p;

void xorshift1024_init64(uint64_t x) {
    int i;
	splitmix64_init(x);
	for(i = 0; i < 16; i++) {
		xorshift1024_s[i] = splitmix64_next();
	}
}

uint64_t xorshift1024_next(void) {
	const uint64_t s0 = xorshift1024_s[xorshift1024_p];
	uint64_t s1 = xorshift1024_s[xorshift1024_p = (xorshift1024_p + 1) & 15];
	s1 ^= s1 << 31; // a
	xorshift1024_s[xorshift1024_p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
	return xorshift1024_s[xorshift1024_p] * UINT64_C(1181783497276652981);
}


/* This is the jump function for the generator. It is equivalent
   to 2^512 calls to next(); it can be used to generate 2^512
   non-overlapping subsequences for parallel computations. */

void xorshift1024_jump(void) {
	static const uint64_t JUMP[] = { 0x84242f96eca9c41d,
		0xa3c65b8776f96855, 0x5b34a39f070b5837, 0x4489affce4f31a1e,
		0x2ffeeb0a48316f40, 0xdc2d9891fe68c022, 0x3659132bb12fea70,
		0xaac17d8efa43cab8, 0xc4cb815590989b13, 0x5ee975283d71c93b,
		0x691548c86c1bd540, 0x7910c41d10a1e6a5, 0x0b5fc64563b3e2a8,
		0x047f7684e9fc949d, 0xb99181f2d8f685ca, 0x284600e3f30e38c3
	};
    int i,j,b;

	uint64_t t[16] = { 0 };
	for(i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b)
				for(j = 0; j < 16; j++)
					t[j] ^= xorshift1024_s[(j + xorshift1024_p) & 15];
			xorshift1024_next();
		}

	for(j = 0; j < 16; j++)
		xorshift1024_s[(j + xorshift1024_p) & 15] = t[j];
}
