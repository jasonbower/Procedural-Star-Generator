#ifndef STAR_STANDARD_H
#define STAR_STANDARD_H

typedef void* STAR;

STAR star_init_default(void);

void star_generate_random(STAR hStar);

void star_print_details(STAR hStar);

void star_destroy(STAR* phStar);

#endif