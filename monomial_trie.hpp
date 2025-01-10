#ifndef MONOMIAL_TRIE_H
#define MONOMIAL_TRIE_H

#include <cstdint>
#include <span>
#include <vector>

#include "kommunopp.hpp"

namespace kommunopp {

template<internal::value_concept V = uint8_t, typename I = uint32_t>
struct monomial_trie {
  std::vector<uint64_t> nodes;
  size_t child_size;

  using match = std::pair<I, size_t>;

  monomial_trie(size_t n)
    : child_size(n) {
    nodes.reserve((child_size + 1) * 1024);
    // insert root
    nodes.insert(nodes.end(), child_size + 1, 0);
  }
  // -----------------------------------------------------------------

  void print() {
    size_t i = 0;
    while(i < nodes.size()) {
      std::cout << "nodes[" << i << "] = " << nodes[i++];
      std::cout << " : [";
      for(size_t j = 1; j < child_size + 1; j++)
        std::cout << nodes[i++] << ", ";
      std::cout << "]\n";
    }
  }
  // -----------------------------------------------------------------

  size_t size() { return nodes.size(); }
  // -----------------------------------------------------------------

  void inline insert(const std::span<const V>& entry, I id) {
    size_t cur_idx = 0;

    for(auto c : entry) {
      if(nodes[cur_idx + c] == 0) {
        nodes[cur_idx + c] = nodes.size();
        nodes.insert(nodes.end(), child_size + 1, 0);
      }
      cur_idx = nodes[cur_idx + c];
    }
    assert(nodes[cur_idx] == 0);
    nodes[cur_idx] = id;
  }
  // -----------------------------------------------------------------

  void inline insert_rev(const std::span<const V>& entry, I id) {
    size_t cur_idx = 0;

    for(auto it = entry.rbegin(); it != entry.rend(); it++) {
      auto c = *it;
      if(nodes[cur_idx + c] == 0) {
        nodes[cur_idx + c] = nodes.size();
        nodes.insert(nodes.end(), child_size + 1, 0);
      }
      cur_idx = nodes[cur_idx + c];
    }
    nodes[cur_idx] = id;
  }
  // -----------------------------------------------------------------
  size_t inline starts_sequence(size_t idx, const std::span<const V>& word) {
    uint64_t cur_idx = idx;
    for(auto c : word) {
      cur_idx = nodes[cur_idx + c];
      if(cur_idx == 0)
        break;
    }
    return cur_idx;
  }
  // -----------------------------------------------------------------

  void inline collect_words(size_t idx,
                            std::vector<match>& res,
                            size_t depth) const {
    uint64_t word = nodes[idx];
    if(word != 0)
      res.emplace_back(word, depth);

    // Recursively visit all children
    for(size_t i = idx + 1; i < idx + 1 + child_size; i++) {
      uint64_t child_idx = nodes[i];
      if(child_idx != 0)
        collect_words(child_idx, res, depth);
    }
  }
  // -----------------------------------------------------------------
  void inline collect_words_adaptive(size_t idx,
                                     std::vector<match>& res,
                                     size_t depth) const {

    uint64_t word = nodes[idx];
    if(word != 0)
      res.emplace_back(word, depth);

    // Recursively visit all children
    for(size_t i = idx + 1; i < idx + 1 + child_size; i++) {
      uint64_t child_idx = nodes[i];
      if(child_idx != 0)
        collect_words_adaptive(child_idx, res, depth + 1);
    }
  }
  // -----------------------------------------------------------------
  void inline inclusions(const std::span<const V>& word,
                         std::vector<match>& inclusions) {

    size_t end_idx;
    for(size_t idx = 0; idx < nodes.size(); idx += child_size + 1) {
      end_idx = starts_sequence(idx, word);
      if(end_idx == 0)
        continue;
      collect_words_adaptive(end_idx, inclusions, 0);
    }
  }

  // -----------------------------------------------------------------

  void inline overlaps_rev(const std::span<const V>& word,
                           std::vector<match>& overlaps) {

    long D = static_cast<long>(word.size());
    for(long d = 1; d < D; d++) {
      uint64_t cur_idx = 0;
      // iterate over suffix of length d
      // but in reversed order
      for(auto it = word.rbegin() + static_cast<long>(d); it != word.rend();
          it++) {
        cur_idx = nodes[cur_idx + *it];
        if(cur_idx == 0)
          break;
      }
      if(cur_idx) {
        // compute overlaps, but skip current node
        // as this is an inclusion
        for(size_t i = 1; i < child_size + 1; i++) {
          uint64_t child_idx = nodes[cur_idx + i];
          if(child_idx != 0) {
            assert(D >= d);
            collect_words(child_idx,
                          overlaps,
                          static_cast<size_t>(D) - static_cast<size_t>(d));
          }
        }
      }
    }
  }
  // -----------------------------------------------------------------

  size_t inline prefixes(const std::span<const V>& prefix,
                         std::vector<match>& res,
                         size_t depth) {
    size_t cur_idx = 0;

    // check if prefix actually appears
    for(auto c : prefix) {
      uint64_t word = nodes[cur_idx];
      if(word)
        res.emplace_back(word, depth);
      cur_idx = nodes[cur_idx + c];
      // no more matching prefixes
      if(cur_idx == 0)
        return cur_idx;
    }
    if(nodes[cur_idx])
      res.emplace_back(nodes[cur_idx], depth);
    return cur_idx;
  }
  // -----------------------------------------------------------------

  inline const std::vector<match>& divisors(
    const std::span<const V>& dividend) {
    size_t D = dividend.size();
    divisors_res.clear();
    divisors_res.reserve(D);
    for(size_t d = 1; d <= D; d++)
      prefixes(dividend.last(d), divisors_res, D - d);
    return divisors_res;
  }
  // -----------------------------------------------------------------
  void inline overlaps_and_inclusions(const std::span<const V>& word,
                                      std::vector<match>& overlaps,
                                      std::vector<match>& inclusions) {

    size_t D = word.size();
    for(size_t d = 1; d <= D; d++) {

      // compute inclusions
      uint64_t cur_idx = prefixes(word.last(d), inclusions, D - d);
      if(cur_idx == 0)
        continue;

      // compute overlaps, but skip current node
      // as this is an inclusion
      for(size_t i = cur_idx + 1; i < cur_idx + 1 + child_size; i++) {
        uint64_t child_idx = nodes[i];
        if(child_idx)
          collect_words(child_idx, overlaps, d);
      }
    }
  }

  private:
  std::vector<match> divisors_res;
};
}

#endif// MONOMIAL_TRIE_H
