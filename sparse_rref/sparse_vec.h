#ifndef SPARSE_VEC_H
#define SPARSE_VEC_H

#include <iostream>
#include "scalar.h"

// Memory management

template <typename T>
inline T* s_malloc(const size_t size) {
  return (T*)std::malloc(size * sizeof(T));
}


template <typename T>
inline void s_free(T* s) {
	std::free(s);
}

template <typename T>
inline T* s_realloc(T* s, const size_t size) {
  return (T*)std::realloc(s, size * sizeof(T));
}

template<typename T>
struct sparse_vec_struct {
  ulong nnz = 0;
  ulong alloc = 0;
  ulong* indices = NULL;
  T* entries = NULL;
};

template<typename T>
using sparse_vec_t = struct sparse_vec_struct<T>[1];

typedef sparse_vec_t<ulong> snmod_vec_t;
typedef sparse_vec_t<fmpz> sfmpz_vec_t;

// sparse_vec

// memory management
template<typename T>
void
sparse_vec_realloc(sparse_vec_t<T> vec, ulong alloc) {
  if(alloc == vec->alloc)
    return;
  // so sparse_vec_realloc(vec,vec->alloc) is useless
  ulong old_alloc = vec->alloc;
  vec->alloc = alloc;
  if(vec->alloc > old_alloc) {
    // enlarge: init later
    vec->indices = s_realloc(vec->indices, vec->alloc);
    vec->entries = s_realloc(vec->entries, vec->alloc);

    if constexpr(std::is_same_v<T, fmpz>) {
      for(ulong i = old_alloc; i < vec->alloc; i++)
        fmpz_init((fmpz*)(vec->entries) + i);
    }
  } else {
    // shrink: clear first
    if constexpr(std::is_same_v<T, fmpz>) {
      for(ulong i = vec->alloc; i < old_alloc; i++)
        fmpz_clear((fmpz*)(vec->entries) + i);
    }
    vec->indices = s_realloc(vec->indices, vec->alloc);
    vec->entries = s_realloc(vec->entries, vec->alloc);
  }
}

#define sparse_vec_entry_pointer(vec, index) ((vec)->entries + (index))

// alloc at least 1 to make sure that indices and entries are not NULL
template<typename T>
inline void
sparse_vec_init(sparse_vec_t<T> vec, ulong alloc = 1) {
  vec->nnz = 0;
  vec->alloc = alloc;
  vec->indices = s_malloc<ulong>(vec->alloc);
  vec->entries = s_malloc<T>(alloc);
  if constexpr(std::is_same_v<T, fmpz>) {
    for(ulong i = 0; i < alloc; i++)
      fmpz_init(vec->entries + i);
  }
}


// set zero and clear memory
template<typename T>
inline void
sparse_vec_clear(sparse_vec_t<T> vec) {
  vec->nnz = 0;
  vec->alloc = 0;
  s_free(vec->indices);
  vec->indices = NULL;
  if constexpr(std::is_same_v<T, fmpz>) {
    for(size_t i = 0; i < vec->alloc; i++)
      fmpz_clear(vec->entries + i);
  }
  s_free(vec->entries);
  vec->entries = NULL;
}

template<typename T>
inline T*
sparse_vec_entry(sparse_vec_t<T> vec, ulong index, const bool isbinary = true) {
  if(vec->nnz == 0 || index < vec->indices[0]
     || index > vec->indices[vec->nnz - 1])
    return NULL;
  ulong* ptr;
  if(isbinary)
    ptr = binarysearch(vec->indices, vec->indices + vec->nnz, index);
  else
    ptr = std::find(vec->indices, vec->indices + vec->nnz, index);
  if(ptr == vec->indices + vec->nnz)
    return NULL;
  return sparse_vec_entry_pointer(vec, ptr - vec->indices);
}

template<typename T>
void
sparse_vec_canonicalize(sparse_vec_t<T> vec) {
  ulong new_nnz = 0;
  ulong i = 0;
  for(; i < vec->nnz; i++) {
    if(!scalar_is_zero(sparse_vec_entry_pointer(vec, i)))
      break;
  }
  for(; i < vec->nnz; i++) {
    if(scalar_is_zero(sparse_vec_entry_pointer(vec, i)))
      continue;
    vec->indices[new_nnz] = vec->indices[i];
    scalar_set(sparse_vec_entry_pointer(vec, new_nnz),
               sparse_vec_entry_pointer(vec, i));
    new_nnz++;
  }
  vec->nnz = new_nnz;
}

template<typename T>
inline void
sparse_vec_compress(sparse_vec_t<T> vec) {
  sparse_vec_canonicalize(vec);
  sparse_vec_realloc(vec, vec->nnz);
}

// debug only, not used to the large vector
template<typename T>
void
print_vec_info(const sparse_vec_t<T> vec) {
  std::cout << "-------------------" << std::endl;
  std::cout << "nnz: " << vec->nnz << std::endl;
  std::cout << "alloc: " << vec->alloc << std::endl;
  std::cout << "indices: ";
  for(size_t i = 0; i < vec->nnz; i++)
    std::cout << vec->indices[i] << " ";
  std::cout << "\nentries: ";
  for(size_t i = 0; i < vec->nnz; i++)
    std::cout << std::to_string(vec->entries + i) << " ";
  std::cout << std::endl;
}

#endif
