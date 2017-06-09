#ifndef SPLITMIX64_H
#define SPLITMIX64_H

void splitmix64_init(uint64_t x);
uint64_t splitmix64_next(void);

#endif
