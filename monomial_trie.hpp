#ifndef MONOMIAL_TRIE_H
#define MONOMIAL_TRIE_H

#include <cstdint>
#include <span>
#include <vector>

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
}

#endif// MONOMIAL_TRIE_H
