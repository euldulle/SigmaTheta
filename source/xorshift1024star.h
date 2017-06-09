#ifndef XORSHIFT1024_H
#define XORSHIFT1024_H

void xorshift1024_init64(uint64_t x);
uint64_t xorshift1024_next(void);
void xorshift1024_jump(void);

#endif
