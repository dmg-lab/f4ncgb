#ifndef SPARSE_VEC_H
#define SPARSE_VEC_H

#include "flint/fmpz.h"
#include <algorithm>
#include <iostream>

// Memory management

template<typename T>
inline T*
s_malloc(const size_t size) {
  return (T*)std::malloc(size * sizeof(T));
}

template<typename T>
inline void
s_free(T* s) {
  std::free(s);
}

template<typename T>
inline T*
s_realloc(T* s, const size_t size) {
  assert(size > 0);
  return (T*)std::realloc(s, size * sizeof(T));
}

// Scalar basics

template<typename T>
inline T*
binarysearch(T* begin, T* end, T val) {
  auto ptr = std::lower_bound(begin, end, val);
  if(ptr == end || *ptr == val)
    return ptr;
  else
    return end;
}

template<typename T>
struct sparse_vec_struct {
  bool is_new_piv = false;
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

  vec->alloc = alloc;
  // enlarge: init later
  vec->indices = s_realloc(vec->indices, vec->alloc);
  vec->entries = s_realloc(vec->entries, vec->alloc);
}

#define sparse_vec_entry_pointer(vec, index) ((vec)->entries + (index))

template<typename T>
inline void
sparse_vec_init(sparse_vec_t<T> vec, ulong alloc = 0) {
  vec->is_new_piv = false;
  vec->nnz = 0;
  vec->alloc = alloc;
  vec->indices = s_malloc<ulong>(alloc);
  vec->entries = s_malloc<T>(alloc);
}

// set zero and clear memory
template<typename T>
inline void
sparse_vec_clear(sparse_vec_t<T> vec) {
  s_free(vec->indices);
  s_free(vec->entries);
  vec->indices = NULL;
  vec->entries = NULL;
  vec->nnz = 0;
  vec->alloc = 0;
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
    std::cout << *(vec->entries + i) << " ";
  std::cout << std::endl;
}

#endif
