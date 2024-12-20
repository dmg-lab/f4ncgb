/*
* Copyright (C) 2018 Christopher Gilbert.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#ifndef AHO_CORASICK_HPP
#define AHO_CORASICK_HPP

#include <algorithm>
#include <cctype>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <queue>
#include <utility>
#include <vector>

namespace aho_corasick {

	// class interval
	class interval {
		size_t d_start;
		size_t d_end;

	public:
		interval(size_t start, size_t end)
			: d_start(start)
			, d_end(end) {}

		size_t get_start() const { return d_start; }
		size_t get_end() const { return d_end; }
		size_t size() const { return d_end - d_start + 1; }

		bool overlaps_with(const interval& other) const {
			return d_start <= other.d_end && d_end >= other.d_start;
		}

		bool overlaps_with(size_t point) const {
			return d_start <= point && point <= d_end;
		}

		bool operator <(const interval& other) const {
			return get_start() < other.get_start();
		}

		bool operator !=(const interval& other) const {
			return get_start() != other.get_start() || get_end() != other.get_end();
		}

		bool operator ==(const interval& other) const {
			return get_start() == other.get_start() && get_end() == other.get_end();
		}
	};


	// class emit
	template<typename CharType>
	class emit: public interval {
	public:
		typedef std::vector<CharType>  string_type;
		typedef std::vector<CharType>& string_ref_type;

	private:
		string_type d_keyword;
		unsigned    d_index = 0;

	public:
		emit()
			: interval(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
			, d_keyword() {}

		emit(size_t start, size_t end, string_type keyword, unsigned index)
			: interval(start, end)
			, d_keyword(keyword), d_index(index) {}

		string_type get_keyword() const { return string_type(d_keyword); }
		unsigned get_index() const { return d_index; }
		bool is_empty() const { return (get_start() == std::numeric_limits<size_t>::max() && get_end() == std::numeric_limits<size_t>::max()); }
	};

	// class token
	template<typename CharType>
	class token {
	public:
		enum token_type{
			TYPE_FRAGMENT,
			TYPE_MATCH,
		};

		using string_type     = std::vector<CharType>;
		using string_ref_type = std::vector<CharType>&;
		using emit_type       = emit<CharType>;

	private:
		token_type  d_type;
		string_type d_fragment;
		emit_type   d_emit;

	public:
		explicit token(string_ref_type fragment)
			: d_type(TYPE_FRAGMENT)
			, d_fragment(fragment)
			, d_emit() {}

		token(string_ref_type fragment, const emit_type& e)
			: d_type(TYPE_MATCH)
			, d_fragment(fragment)
			, d_emit(e) {}

		bool is_match() const { return (d_type == TYPE_MATCH); }
		string_type get_fragment() const { return string_type(d_fragment); }
		emit_type get_emit() const { return d_emit; }
	};

	// class state
	template<typename CharType>
	class state {
	public:
		typedef state<CharType>*                 ptr;
		typedef std::unique_ptr<state<CharType>> unique_ptr;
		typedef std::vector<CharType>      string_type;
		typedef std::vector<CharType>&     string_ref_type;
		typedef std::pair<string_type, unsigned> key_index;
		typedef std::set<key_index>              string_collection;
		typedef std::vector<ptr>                 state_collection;
		typedef std::vector<CharType>            transition_collection;

	private:
		size_t                         d_depth;
		ptr                            d_root;
		std::map<CharType, unique_ptr> d_success;
		ptr                            d_failure;
		string_collection              d_emits;
		std::set<string_type>          keywords;

	public:
		state(): state(0) {}

		explicit state(size_t depth)
			: d_depth(depth)
			, d_root(depth == 0 ? this : nullptr)
			, d_success()
			, d_failure(nullptr)
			, d_emits() {}

		ptr next_state(CharType character) const {
			return next_state(character, false);
		}

		ptr next_state_ignore_root_state(CharType character) const {
			return next_state(character, true);
		}

		ptr add_state(CharType character) {
			auto next = next_state_ignore_root_state(character);
			if (next == nullptr) {
				next = new state<CharType>(d_depth + 1);
				d_success[character].reset(next);
			}
			return next;
		}

		size_t get_depth() const { return d_depth; }

		void add_emit(string_ref_type keyword, unsigned index) {
			d_emits.insert(std::make_pair(keyword, index));
		}

		void add_emit(const string_collection& emits) {
			for (const auto& e : emits) {
				string_type str(e.first);
				add_emit(str, e.second);
			}
		}
		
		void add_keyword(string_ref_type key) {
		    keywords.insert(key);
		}
		
		std::set<string_type> get_keywords() const { return keywords; }

		string_collection get_emits() const { return d_emits; }

		ptr failure() const { return d_failure; }

		void set_failure(ptr fail_state) { d_failure = fail_state; }

		state_collection get_states() const {
			state_collection result;
			for (auto it = d_success.cbegin(); it != d_success.cend(); ++it) {
				result.push_back(it->second.get());
			}
			return state_collection(result);
		}

		transition_collection get_transitions() const {
			transition_collection result;
			for (auto it = d_success.cbegin(); it != d_success.cend(); ++it) {
				result.push_back(it->first);
			}
			return transition_collection(result);
		}

	private:
		ptr next_state(CharType character, bool ignore_root_state) const {
			ptr result = nullptr;
			auto found = d_success.find(character);
			if (found != d_success.end()) {
				result = found->second.get();
			} else if (!ignore_root_state && d_root != nullptr) {
				result = d_root;
			}
			return result;
		}
	};

	template<typename CharType>
	class trie {
	public:
		using string_type = std::vector < CharType > ;
		using string_ref_type = std::vector<CharType>&;

		typedef state<CharType>         state_type;
		typedef state<CharType>*        state_ptr_type;
		typedef token<CharType>         token_type;
		typedef emit<CharType>          emit_type;
		typedef std::vector<token_type> token_collection;
		typedef std::vector<emit_type>  emit_collection;

		class config {
			bool d_allow_overlaps;
			bool d_only_whole_words;
			bool d_case_insensitive;

		public:
			config()
				: d_allow_overlaps(true)
				, d_only_whole_words(false)
				, d_case_insensitive(false) {}
		};

	private:
		std::unique_ptr<state_type> d_root;
		config                      d_config;
		bool                        d_constructed_failure_states;
		unsigned                    d_num_keywords = 0;

	public:
		trie(): trie(config()) {}

		trie(const config& c)
			: d_root(new state_type())
			, d_config(c)
			, d_constructed_failure_states(false) {}
		
		void insert(string_type keyword, bool reverse=false) {
			if (keyword.empty())
				return;
			state_ptr_type cur_state = d_root.get();
			for (const auto& ch : keyword) {
				cur_state = cur_state->add_state(ch);
			}
			cur_state->add_emit(keyword, d_num_keywords++);
			if (reverse)
			    std::reverse(keyword.begin(), keyword.end());
			cur_state->add_keyword(keyword);
			d_constructed_failure_states = false;
		}

		template<class InputIterator>
		void insert(InputIterator first, InputIterator last) {
			for (InputIterator it = first; first != last; ++it) {
				insert(*it);
			}
		}
		
		emit_collection parse_text(const std::span<const CharType>::iterator& begin, 
		                           const std::span<const CharType>::iterator& end) {
			check_construct_failure_states();
			size_t pos = 0;
			state_ptr_type cur_state = d_root.get();
			emit_collection collected_emits;
			for (auto it = begin; it != end; it++) {
				cur_state = get_state(cur_state, *it);
				store_emits(pos, cur_state, collected_emits);
				pos++;
			}
			return collected_emits;
		}
		
		template <typename Iterator>
		emit_collection prefixes(const Iterator& begin, 
		                         const Iterator& end) {
			emit_collection collected_emits;
						
			// check if prefix is in trie
			state_ptr_type cur_state = d_root.get();
			cur_state = starts_sequence(cur_state, begin, end);

			// traverse remaining trie and store all keywords
			if (cur_state != nullptr)
                all_keywords_in_subtrie(cur_state, collected_emits);
            return collected_emits;
		}
				
        emit_collection inclusions(const std::span<const CharType>::iterator& begin, 
		                         const std::span<const CharType>::iterator& end) {
			emit_collection collected_emits;
						
			state_ptr_type cur_state = d_root.get();
            std::queue<state_ptr_type> q;
            q.push(cur_state);

			while (!q.empty()) {
				cur_state = q.front();
				cur_state = starts_sequence(cur_state, begin, end);
				
				if (cur_state != nullptr)
				    all_keywords_in_subtrie(cur_state, collected_emits);

                // append all next states
                cur_state = q.front();
				for (const auto& t : cur_state->get_transitions())
				    q.push(cur_state->next_state(t));
				q.pop();
			}			
			return collected_emits;
		}



	private:
	    template <typename Iterator>
	    inline state_ptr_type starts_sequence(state_ptr_type & cur_state,
	                                        const Iterator& begin, 
		                                    const Iterator& end) {
	        for (auto it = begin; it != end and cur_state != nullptr; it++)
			    cur_state = cur_state->next_state_ignore_root_state(*it);
			return cur_state;

	    }
	    
	    inline void all_keywords_in_subtrie(state_ptr_type & cur_state, emit_collection & collected_emits) {
		    std::queue<std::pair<state_ptr_type,size_t>> q;
            q.push( std::make_pair(cur_state,0) );

			while (!q.empty()) {
				auto [cur_state,pos] = q.front();
				store_keywords(pos, cur_state, collected_emits);
				for (const auto& t : cur_state->get_transitions()) {
				    q.push( std::make_pair(cur_state->next_state(t),pos+1) );
				}
				q.pop();
			}			
		}

		state_ptr_type get_state(state_ptr_type cur_state, CharType c) const {
			state_ptr_type result = cur_state->next_state(c);
			while (result == nullptr) {
				cur_state = cur_state->failure();
				result = cur_state->next_state(c);
			}
			return result;
		}

		void check_construct_failure_states() {
			if (!d_constructed_failure_states) {
				construct_failure_states();
			}
		}

		void construct_failure_states() {
			std::queue<state_ptr_type> q;
			for (auto& depth_one_state : d_root->get_states()) {
				depth_one_state->set_failure(d_root.get());
				q.push(depth_one_state);
			}
			d_constructed_failure_states = true;

			while (!q.empty()) {
				auto cur_state = q.front();
				for (const auto& transition : cur_state->get_transitions()) {
					state_ptr_type target_state = cur_state->next_state(transition);
					q.push(target_state);

					state_ptr_type trace_failure_state = cur_state->failure();
					while (trace_failure_state->next_state(transition) == nullptr) {
						trace_failure_state = trace_failure_state->failure();
					}
					state_ptr_type new_failure_state = trace_failure_state->next_state(transition);
					target_state->set_failure(new_failure_state);
					target_state->add_emit(new_failure_state->get_emits());
				}
				q.pop();
			}
		}

		void store_emits(size_t pos, state_ptr_type cur_state, emit_collection& collected_emits) const {
			auto emits = cur_state->get_emits();
			for (const auto& str : emits) {
				auto emit_str = typename emit_type::string_type(str.first);
				collected_emits.push_back(emit_type(pos - emit_str.size() + 1, pos, emit_str, str.second));
			}
		}
		
		void store_keywords(size_t pos, state_ptr_type cur_state, emit_collection& collected_emits) const {
			auto keywords = cur_state->get_keywords();
			for (const auto& str : keywords) {
			    auto emit_str = typename emit_type::string_type(str);
				collected_emits.push_back(emit_type(pos - emit_str.size() + 1, pos, emit_str, 0));
			}
		}
	};

} // namespace aho_corasick

#endif // AHO_CORASICK_HPP