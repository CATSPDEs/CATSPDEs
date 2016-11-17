#pragma once
#include <list>

// ssize_t
#if __BITS_PER_LONG != 64
	typedef int ssize_t;
#else
	typedef long ssize_t;
#endif

typedef  size_t Index;
typedef ssize_t SignedIndex;
typedef unsigned short LocalIndex;

typedef std::list<Index>      Indicies;
typedef std::list<LocalIndex> LocalIndicies;
