#include "monomial_trie.hpp"

namespace kommunopp {
// We instantiate this variant directly, so that the other compilation unit does
// not have to reproduce the same class tens of times.
template struct kommunopp::monomial_trie<uint8_t, uint32_t>;
}
