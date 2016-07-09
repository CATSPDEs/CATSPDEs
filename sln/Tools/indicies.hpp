#pragma once
#include <list>

using namespace std;

// ssize_t
#if __BITS_PER_LONG != 64
typedef int ssize_t;
#else
typedef long ssize_t;
#endif

typedef  size_t Index;
typedef ssize_t SignedIndex;
typedef unsigned short LocalIndex;

typedef list<Index>      Indicies;
typedef list<LocalIndex> LocalIndicies;
