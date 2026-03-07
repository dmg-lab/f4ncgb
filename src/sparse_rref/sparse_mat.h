#ifndef SPARSE_MAT_H
#define SPARSE_MAT_H

#include "sparse_vec.h"
#include <boost/multiprecision/gmp.hpp>

template<typename T>
struct sparse_mat_struct {
  ulong nrow;
  ulong ncol;
  sparse_vec_struct<T>* rows;
};

template<typename T>
using sparse_mat_t = struct sparse_mat_struct<T>[1];

typedef sparse_mat_t<ulong> snmod_mat_t;
typedef sparse_mat_t<fmpz> sfmpz_mat_t;

#define sparse_mat_row(mat, ind) ((mat)->rows + (ind))

template<typename T>
inline void
sparse_mat_init(sparse_mat_t<T> mat,
                ulong nrow,
                ulong ncol) {
  mat->nrow = nrow;
  mat->ncol = ncol;
  mat->rows = s_malloc<sparse_vec_struct<T>>(nrow);
}

template<typename T>
inline void
sparse_mat_clear(sparse_mat_t<T> mat) {  
  for(size_t i = 0; i < mat->nrow; i++)
    sparse_vec_clear(sparse_mat_row(mat, i));
  s_free(mat->rows);
  mat->nrow = 0;
  mat->ncol = 0;
  mat->rows = NULL;
}

template<typename T>
inline sparse_vec_struct<T>*
sparse_mat_row_init(sparse_mat_t<T> mat, ulong i, ulong alloc) {
  sparse_vec_init(sparse_mat_row(mat, i), alloc);
  return sparse_mat_row(mat, i);
}

template<typename T>
inline ulong
sparse_mat_nnz(const sparse_mat_t<T> mat) {
  ulong nnz = 0;
  for(size_t i = 0; i < mat->nrow; i++)
    nnz += sparse_mat_row(mat, i)->nnz;
  return nnz;
}

template<typename T>
inline T*
sparse_mat_entry(sparse_mat_t<T> mat,
                 ulong row,
                 ulong col,
                 bool isbinary = true) {
  return sparse_vec_entry(sparse_mat_row(mat, row), col, isbinary);
}

template<typename T, typename S>
void
sparse_mat_write(sparse_mat_t<T> mat, S& st) {
  st << mat->nrow << ' ' << mat->ncol << ' ' << sparse_mat_nnz(mat) << '\n';
  for(size_t i = 0; i < mat->nrow; i++) {
    for(size_t j = 0; j < mat->ncol; j++) {
      auto c = sparse_mat_entry(mat, i, j);
      boost::multiprecision::mpz_int cc;
      if(c == nullptr)
        st << 0 << ' ';
      else if constexpr(std::is_same_v<T, fmpz>) {
        fmpz_get_mpz(cc.backend().data(), c);
        st << cc << ' ';
      } else
        st << (int)(*c) << ' ';
    }
    st << "\n";
  }
}

#endif
