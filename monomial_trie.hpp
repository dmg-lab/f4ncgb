#ifndef MONOMIAL_TRIE_H
#define MONOMIAL_TRIE_H

#include <algorithm>
#include <chrono>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <unordered_set>
#include <vector>

#include "./aho_corasick/src/aho_corasick/aho_corasick.hpp"

#include "ambiguity.hpp"
#include "kommunopp.hpp"


namespace kommunopp {

template<size_t N, typename I = uint32_t>
struct trienode {
  uint64_t children[N + 1] = {};
  I is_word = 0;
};

template<size_t N, internal::value_concept V = uint8_t, typename I = uint32_t>
struct monomial_trie {
  std::vector<trienode<N, I>> nodes;
  size_t child_size = N + 1;

  using match = std::pair<I, size_t>;

  monomial_trie() {
    nodes.reserve(1024);
    nodes.emplace_back();
  }
  // -----------------------------------------------------------------

  void print() {
    size_t i = 0;
    for(auto n : nodes) {
      std::cout << "nodes[" << i++ << "] = " << (int)n.is_word;
      std::cout << " : [";
      for(auto c : n.children)
        std::cout << (int)c << ", ";
      std::cout << "]\n";
    }
  }
  // -----------------------------------------------------------------

  size_t size() { return nodes.size(); }
  // -----------------------------------------------------------------

  void inline insert(const std::span<const V>& entry, I id) {
    size_t cur_idx = 0;

    for(auto c : entry) {
      if(nodes[cur_idx].children[c] == 0) {
        nodes[cur_idx].children[c] = nodes.size();
        nodes.emplace_back();
      }
      cur_idx = nodes[cur_idx].children[c];
    }
    assert(nodes[cur_idx].is_word == 0);
    nodes[cur_idx].is_word = id;
  }
  // -----------------------------------------------------------------

  void inline insert_rev(const std::span<const V>& entry, I id) {
    size_t cur_idx = 0;

    for(auto it = entry.rbegin(); it != entry.rend(); it++) {
      auto c = *it;
      if(!nodes[cur_idx].children[c]) {
        nodes[cur_idx].children[c] = nodes.size();
        nodes.emplace_back();
      }
      cur_idx = nodes[cur_idx].children[c];
    }
    nodes[cur_idx].is_word = id;
  }
  // -----------------------------------------------------------------
  size_t inline starts_sequence(size_t idx, const std::span<const V>& word) {
    size_t cur_idx = idx;
    for(auto c : word) {
      cur_idx = nodes[cur_idx].children[c];
      if(cur_idx == 0)
        break;
    }
    return cur_idx;
  }
  // -----------------------------------------------------------------

  void inline collect_words(size_t idx,
                            std::vector<match>& res,
                            size_t depth) const {

    auto node = nodes[idx];
    if(node.is_word)
      res.emplace_back(node.is_word, depth);

    // Recursively visit all children
    for(size_t i = 1; i < child_size; i++) {
      I child_idx = node.children[i];
      if(child_idx)
        collect_words(child_idx, res, depth);
    }
  }
  // -----------------------------------------------------------------
  void inline collect_words_adaptive(size_t idx,
                                     std::vector<match>& res,
                                     size_t depth) const {

    auto node = nodes[idx];
    if(node.is_word)
      res.emplace_back(node.is_word, depth);

    // Recursively visit all children
    for(auto c : node.children) {
      if(c)
        collect_words_adaptive(c, res, depth + 1);
    }
  }
  // -----------------------------------------------------------------
  void inline inclusions(const std::span<const V>& word,
                         std::vector<match>& inclusions) {

    size_t end_idx;
    for(size_t idx = 0; idx < nodes.size(); idx++) {
      end_idx = starts_sequence(idx, word);
      if(end_idx == 0)
        continue;
      collect_words_adaptive(end_idx, inclusions, 0);
    }
  }

  // -----------------------------------------------------------------

  void inline overlaps_rev(const std::span<const V>& word,
                           std::vector<match>& overlaps) {

    size_t D = word.size();
    for(size_t d = 1; d < D; d++) {
      size_t cur_idx = 0;
      // iterate over suffix of length d
      // but in reversed order
      for(auto it = word.rbegin() + d; it != word.rend(); it++) {
        cur_idx = nodes[cur_idx].children[*it];
        if(cur_idx == 0)
          break;
      }
      if(cur_idx) {
        // compute overlaps, but skip current node
        // as this is an inclusion
        for(auto c : nodes[cur_idx].children) {
          if(c)
            collect_words(c, overlaps, D - d);
        }
      }
    }
  }
  // -----------------------------------------------------------------

  size_t prefixes(const std::span<const V>& prefix,
                  std::vector<match>& res,
                  size_t depth) {
    size_t cur_idx = 0;

    // check if prefix actually appears
    for(auto c : prefix) {
      auto n = nodes[cur_idx];
      if(n.is_word)
        res.emplace_back(n.is_word, depth);
      cur_idx = n.children[c];
      // no more matching prefixes
      if(cur_idx == 0)
        return cur_idx;
    }
    if(nodes[cur_idx].is_word)
      res.emplace_back(nodes[cur_idx].is_word, depth);
    return cur_idx;
  }
  // -----------------------------------------------------------------

  std::vector<match> inline divisors(const std::span<const V>& dividend) {
    std::vector<match> res;
    size_t D = dividend.size();
    res.reserve(D);
    for(size_t d = 1; d <= D; d++)
      prefixes(dividend.last(d), res, D - d);
    return res;
  }
  // -----------------------------------------------------------------
  void inline overlaps_and_inclusions(const std::span<const V>& word,
                                      std::vector<match>& overlaps,
                                      std::vector<match>& inclusions) {

    size_t D = word.size();
    for(size_t d = 1; d <= D; d++) {

      // compute inclusions
      auto cur_idx = prefixes(word.last(d), inclusions, D - d);
      if(cur_idx == 0)
        continue;

      // compute overlaps, but skip current node
      // as this is an inclusion
      for(auto c : nodes[cur_idx].children) {
        if(c)
          collect_words(c, overlaps, d);
      }
    }
  }
};

/* template<internal::metadata_concept M = internal::metadata<uint8_t>, */
/*          internal::value_concept V = uint8_t, */
/*          typename I = uint32_t> */
/* class monomial_trie { */
/*   using monomial_store = internal::monomial_store<M, V, I>; */
/*   using ambiguity = ambiguity<I>; */
/*   using amb_hash = ambiguity_hash<I>; */
/*   using monomial = std::span<const V>; */

/*   public: */
/*   aho_corasick::trie<V> T; */

/*   monomial_trie() {} */
/*   ~monomial_trie() = default; */

/*   inline void insert_monomial(const monomial& m) { */
/*     std::vector<V> vec(m.begin(), m.end()); */
/*     T.insert(vec); */
/*   } */
/*   // ----------------------------------------------------------------- */
/*   inline void insert_monomial_reversed(const monomial& m) { */
/*     std::vector<V> vec(m.rbegin(), m.rend()); */
/*     T.insert(vec, true); */
/*   } */
/*   // ----------------------------------------------------------------- */
/*   std::unordered_set<I> compute_divisors(monomial_store& s, I a) { */
/*     monomial m = s[a]; */
/*     auto results = T.parse_text(m.begin(), m.end()); */

/*     std::unordered_set<I> divisors; */
/*     divisors.reserve(results.size()); */

/*     for(auto& r : results) { */
/*       auto vec = r.get_keyword(); */
/*       divisors.insert(s.getid(vec)); */
/*     } */
/*     return divisors; */
/*   } */
/*   // ----------------------------------------------------------------- */
/*   std::unordered_set<ambiguity, amb_hash> compute_overlaps(monomial_store& s, */
/*                                                            I i) { */
/*     /\* */
/*      * Find all overlaps ABC with a = AB */
/*      *\/ */
/*     monomial ab = s[i]; */
/*     std::unordered_set<ambiguity, amb_hash> overlaps; */

/*     for(I id = 1; id < ab.size(); id++) { */
/*       auto results = T.prefixes(ab.begin() + id, ab.end()); */

/*       // compute the ambiguities */
/*       for(auto& r : results) { */
/*         std::vector<V> v = r.get_keyword(); */
/*         monomial bc(v.begin(), v.size()); */

/*         I d = id + v.size(); */
/*         I j = s.getid(bc); */
/*         I aj = s.getid(ab.first(id)); */
/*         I ci = s.getid(bc.last(r.get_end())); */

/*         ambiguity a(d, i, j, 0, ci, aj, 0); */
/*         overlaps.insert(a); */
/*       } */
/*     } */
/*     return overlaps; */
/*   } */
/*   // ----------------------------------------------------------------- */
/*   std::unordered_set<ambiguity, amb_hash> compute_overlaps_reversed( */
/*     monomial_store& s, */
/*     I i) { */
/*     /\* */
/*      * Find all overlaps ABC with a = BC */
/*      *\/ */
/*     monomial bc = s[i]; */
/*     std::unordered_set<ambiguity, amb_hash> overlaps; */

/*     for(I id = 1; id < bc.size(); id++) { */
/*       auto results = T.prefixes(bc.rbegin() + id, bc.rend()); */

/*       // compute the ambiguities */
/*       for(auto& r : results) { */
/*         std::vector<V> v = r.get_keyword(); */
/*         monomial ab(v.begin(), v.size()); */

/*         I d = id + v.size(); */
/*         I j = s.getid(ab); */
/*         I ai = s.getid(ab.first(r.get_end())); */
/*         I cj = s.getid(bc.last(id)); */

/*         ambiguity a(d, j, i, 0, cj, ai, 0); */
/*         overlaps.insert(a); */
/*       } */
/*     } */
/*     return overlaps; */
/*   } */
/*   // ----------------------------------------------------------------- */
/*   std::unordered_set<ambiguity, amb_hash> compute_inclusions(monomial_store& s, */
/*                                                              I i) { */
/*     monomial m = s[i]; */

/*     // Find all inclusions ABC with m = ABC */
/*     std::unordered_set<ambiguity, amb_hash> inclusions; */
/*     auto results = T.parse_text(m.begin(), m.end()); */

/*     for(const auto& r : results) { */
/*       std::vector<V> v = r.get_keyword(); */
/*       monomial b(v.begin(), v.size()); */
/*       // pos say where A ends in ABC */
/*       I pos = r.get_start(); */
/*       I d = m.size(); */
/*       I j = s.getid(b); */

/*       if(i == j) */
/*         continue; */

/*       I aj = s.getid(m.first(pos)); */
/*       I cj = s.getid(m.last(d - pos - b.size())); */

/*       ambiguity a(d, j, i, aj, cj, 0, 0); */
/*       inclusions.insert(a); */
/*     } */

/*     // Find all inclusions ABC with m = B */
/*     auto start = std::chrono::high_resolution_clock().now(); */
/*     results = T.inclusions(m.begin(), m.end()); */
/*     auto end = std::chrono::high_resolution_clock().now(); */
/*     std::chrono::duration<double> elapsed = end - start; */
/*     tt += elapsed.count(); */

/*     for(const auto& r : results) { */
/*       std::vector<V> v = r.get_keyword(); */
/*       monomial abc(v.begin(), v.size()); */
/*       // pos say how long C is in ABC */
/*       I pos = r.get_end(); */
/*       I d = abc.size(); */
/*       I j = s.getid(abc); */

/*       if(i == j) */
/*         continue; */

/*       I ai = s.getid(abc.first(d - pos - m.size())); */
/*       I ci = s.getid(abc.last(pos)); */

/*       ambiguity a(d, i, j, ai, ci, 0, 0); */
/*       inclusions.insert(a); */
/*     } */

/*     return inclusions; */
/*   } */
/* }; */
}

#endif// MONOMIAL_TRIE_H
