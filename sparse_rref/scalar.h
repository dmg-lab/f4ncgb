#ifndef SCALAR_H
#define SCALAR_H

#include <algorithm>
#include "flint/fmpz.h"

template <typename T> inline T* binarysearch(T* begin, T* end, T val) {
	auto ptr = std::lower_bound(begin, end, val);
	if (ptr == end || *ptr == val)
		return ptr;
	else
		return end;
}

// scalar
static inline bool scalar_is_zero(const fmpz_t a) { return fmpz_is_zero(a); }
static inline bool scalar_is_zero(const uint32_t* a) { return (*a) == 0; }

static inline void scalar_set(fmpz_t a, const fmpz_t b) { fmpz_set(a, b); }
static inline void scalar_set(uint32_t* a, const uint32_t* b) { *a = *b; }


#endif
