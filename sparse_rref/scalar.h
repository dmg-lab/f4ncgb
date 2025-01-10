#ifndef SCALAR_H
#define SCALAR_H

#include <type_traits>

#include "flint/fmpq.h"
#include "flint/nmod.h"
#include "flint/ulong_extras.h"
#include "gmp.h"
#include "sparse_base.h"

template <typename T> inline T* binarysearch(T* begin, T* end, T val) {
	auto ptr = std::lower_bound(begin, end, val);
	if (ptr == end || *ptr == val)
		return ptr;
	else
		return end;
}

// scalar
static inline void scalar_init(fmpz_t a) { fmpz_init(a); }

static inline void scalar_clear(fmpz_t a) { fmpz_clear(a); }

template<typename T>
static inline std::string scalar_to_str(T* a) { return std::to_string(*a); }

static inline bool scalar_is_zero(const fmpz_t a) { return fmpz_is_zero(a); }
static inline bool scalar_is_zero(const uint32_t* a) { return (*a) == 0; }

static inline void scalar_set(fmpz_t a, const fmpz_t b) { fmpz_set(a, b); }
static inline void scalar_set(uint32_t* a, const uint32_t* b) { *a = *b; }
template <typename T>
inline void scalar_set(T* a, const T* b, const ulong rank) {
	for (ulong i = 0; i < rank; i++)
		scalar_set(a + i, b + i);
}


#endif
