#ifndef _TIME_RANDOM_H
#define _TIME_RANDOM_H

#include <limits>
#include <time.h>
#include <stdlib.h>

#undef STDLIB_RANDOM
/* #define STDLIB_RANDOM */

#ifdef STDLIB_RANDOM
#include <stdlib.h>
unsigned Random()
{
	return rand();
}

void SRandom(unsigned Seed)
{
	srand(Seed);
}

#else


#define PRANDMAX INT_MAX

static int a = 0, b = 24, arr[55], initialized = 0;
//unsigned Random(void);
//void SRandom(unsigned Seed);
double GetTime(void);
/*unsigned Random()
{
int t;

if (!initialized)
SRandom(7913);
if (a-- == 0)
a = 54;
if (b-- == 0)
b = 54;
if ((t = arr[a] - arr[b]) < 0)
t += PRANDMAX;
return (arr[a] = t);
}*/

/*void SRandom(unsigned Seed)
{
int i, ii, last, next;

Seed %= PRANDMAX;
arr[0] = last = Seed;
for (next = i = 1; i < 55; i++) {
ii = (21 * i) % 55;
arr[ii] = next;
if ((next = last - next) < 0)
next += PRANDMAX;
last = arr[ii];
}
initialized = 1;
a = 0;
b = 24;
for (i = 0; i < 165; i++)
Random();
}*/

#endif
#define ByteLen 8
int Random(int a, int b)
{
	return (rand() % (b - a + 1)) + a;
}
int* RandomArray(int Dimension)
{
	int *output = (int*)malloc(sizeof(int) * Dimension);
	for(int i = 0; i < Dimension; i++) output[i] = i;
	for(int i = 0; i < Dimension; i++)
	{
		int index = Random(0, Dimension - 1);
		int temp = output[i];
		output[i] = output[index];
		output[index] = temp;
	}
	return output;
}
double GetTime()
{
	return (double) clock() / CLOCKS_PER_SEC;
}
double Random01()
{
	return (double)rand() / (double)(RAND_MAX + 1) ;
}
double Random(double min, double max)
{
	return min + (Random01() * (max - min));
}
#endif