/*
 * Copyright (c) 2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <absl/container/btree_map.h>
#include <absl/container/btree_set.h>
#include <absl/container/flat_hash_set.h>
#include <algorithm>
#include <args.hxx>
#include <array>
#include <atomic>
#include <bit>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <iostream>
#include <iterator>
#include <libbio/accumulator.hh>
#include <libbio/algorithm.hh>
#include <libbio/assert.hh>
#include <libbio/copyable_atomic.hh>
#include <libbio/file_handling.hh>
#include <libbio/fmap.hh>
#include <libbio/join_iterator.hh>
#include <libbio/tuple.hh>
#include <libbio/utility.hh>
#include <limits>
#include <mutex>
#include <numeric>
#include <range/v3/view/filter.hpp>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "args.hpp"
#include "configuration.hpp"
#include "gapmer.hpp"
#include "huddinge_distance.hpp"
#include "math.hpp"
#include "pack_characters.hpp"
#include "packed_read.hpp"
#include "reader_adapter.hpp"
#include "version.hpp"

#if defined(__clang__)
#	define SF_DO_NOT_OPTIMISE __attribute__((optnone))
#endif

#define SF_FORMAT_ITALIC "\033[3m"
#define SF_FORMAT_RESET "\033[0m"

namespace lb = libbio;
namespace rsv = ranges::views;


namespace
{
	//// Helpers.

	template <typename t_tuple, typename... t_args>
	t_tuple make_tuple_with_args(t_args && ... args)
	{
		// Use the null pointer to get the type parameters of the tuple.
		return []<typename... t_types>(std::tuple <t_types...> *, t_args... args_){
			return std::tuple <t_types...>{t_types(std::forward <t_args>(args_)...)...};
		}(static_cast <t_tuple *>(nullptr), std::forward <t_args>(args)...);
	}


	template <std::unsigned_integral t_lhs, std::unsigned_integral t_rhs>
	constexpr inline std::common_type_t <t_lhs, t_rhs> ceil_div(t_lhs n, t_rhs d)
	{
		return n / d + (n % d != 0);
	}


	struct access_map_key
	{
		template <typename t_value>
		decltype(auto) operator()(t_value &&value) const { return value.first; }
	};


	struct access_map_value
	{
		template <typename t_value>
		decltype(auto) operator()(t_value &&value) const { return value.second; }
	};


	// Report the first iteration to the callback.
	struct first_iteration_tag{ constexpr static bool is_first() { return true; } };
	struct subsequent_iteration_tag{ constexpr static bool is_first() { return false; } };

	void for_each_(auto it, auto const end, auto &&cb)
	{
		if (it != end)
		{
			cb(it, first_iteration_tag{});
			++it;
			while (it != end)
			{
				cb(it, subsequent_iteration_tag{});
				++it;
			}
		}
	}


	template <typename t_value>
	t_value update_maximum(std::atomic <t_value> &dst, t_value const value)
	{
		auto prev_value{dst.load(std::memory_order_acquire)};
		while (true)
		{
			// compare_exchange_weak() updates prev_value.
			if (dst.compare_exchange_weak(prev_value, value, std::memory_order_acq_rel))
				return value;

			if (value <= prev_value)
				return prev_value;
		}
	}


	//// Correlations-specific.

	constexpr static bool middle_gap_only{true};
	constexpr static std::uint8_t max_gap{15};
	typedef sf::gapmer <middle_gap_only, max_gap> gapmer_type;
	typedef std::vector <gapmer_type> gapmer_vector;


	typedef std::vector <std::uint64_t> packed_word_vector;


	// Colour type for count_pair. (See below.)
	// Grey and red are not mutually exclusive.
	enum class colour_type : std::uint8_t
	{
		white = 0x0,
		grey = 0x1,
		black = 0x3
	};


	// Since we are short on memory, we store the node colours needed for topological
	// sorting in the most significant bits of the count values.
	class count_pair
	{
	private:
		constexpr static std::uint32_t colour_mask{UINT32_C(0x8000'0000)};
		constexpr static std::uint32_t count_mask{UINT32_C(0x7FFF'FFFF)};

		static_assert(0 == ~(colour_mask ^ count_mask));

	private:
		// We need atomicity in processing topologically sorted nodes.
		// (This step does not involve copying.)
		lb::copyable_atomic <std::uint32_t> m_signal_count{};
		lb::copyable_atomic <std::uint32_t> m_background_count{};

	public:
		constexpr static std::uint64_t max_count() { return count_mask; }

		constexpr count_pair() = default;

		constexpr count_pair(std::uint32_t signal_count, std::uint32_t background_count):
			m_signal_count{signal_count},
			m_background_count{background_count}
		{
		}

		std::int64_t signal_count() const { return m_signal_count.load(std::memory_order_relaxed) & count_mask; }
		std::int64_t background_count() const { return m_background_count.load(std::memory_order_relaxed) & count_mask; }

		double fold_change_(double signal_size, double background_size) const;
		double fold_change(double signal_size, double background_size) const;

		// Safe without applying the mask.
		void increment_signal() { m_signal_count.fetch_add(1, std::memory_order_relaxed); }
		void increment_background() { m_background_count.fetch_add(1, std::memory_order_relaxed); }

		inline colour_type colour() const;
		inline void colour_white();
		inline void set_colour(colour_type cc);

		std::int64_t diff() const { return signal_count() - background_count(); }
		std::int64_t sum() const { return signal_count() + background_count(); }

		count_pair with_pseudocounts() const;

		// Atomic operations.
		inline bool add_grey_colour(std::memory_order order);
		inline bool clear_grey_colour(std::memory_order order);
		inline bool has_grey_colour(std::memory_order order) const;
	};


	static_assert(8 == sizeof(count_pair));


	auto count_pair::with_pseudocounts() const -> count_pair
	{
		return {
			1 + m_signal_count.load(std::memory_order_relaxed),
			1 + m_background_count.load(std::memory_order_relaxed)
		};
	}


	double count_pair::fold_change_(double signal_size, double background_size) const
	{
		return (signal_count() * background_size) / (background_count() * signal_size);
	}


	double count_pair::fold_change(double signal_size, double background_size) const
	{
		return with_pseudocounts().fold_change_(signal_size, background_size);
	}


	// Non-atomic.
	colour_type count_pair::colour() const
	{
		std::uint8_t retval{};
		retval |= (m_signal_count.load(std::memory_order_relaxed) & colour_mask) >> 30;
		retval |= (m_background_count.load(std::memory_order_relaxed) & colour_mask) >> 31;
		return static_cast <colour_type>(retval);
	}


	// Non-atomic.
	void count_pair::colour_white()
	{
		static_assert(lb::to_underlying(colour_type::white) == 0);

		m_signal_count.fetch_and(~colour_mask, std::memory_order_relaxed);
		m_background_count.fetch_and(~colour_mask, std::memory_order_relaxed);
	}


	// Non-atomic.
	void count_pair::set_colour(colour_type cc)
	{
		colour_white();

		auto const cc_{lb::to_underlying(cc)};
		m_signal_count.fetch_or((cc_ & 0x2) << 30, std::memory_order_relaxed);
		m_background_count.fetch_or((cc_ & 0x1) << 31, std::memory_order_relaxed);
	}


	// Atomic.
	bool count_pair::add_grey_colour(std::memory_order order)
	{
		static_assert(lb::to_underlying(colour_type::grey) == 0x1);

		// m_background_count has the lower bit.
		auto const res{m_background_count.fetch_or(std::uint32_t{0x1} << 31U, order)};
		return (not bool(0x1U & std::rotl(res, 1)));
	}


	bool count_pair::clear_grey_colour(std::memory_order order)
	{
		static_assert(lb::to_underlying(colour_type::grey) == 0x1);

		// m_background_count has the lower bit.
		auto const res{m_background_count.fetch_and(count_mask, order)};
		return (not bool(0x1U & std::rotl(res, 1)));
	}


	// Atomic.
	bool count_pair::has_grey_colour(std::memory_order order) const
	{
		static_assert(lb::to_underlying(colour_type::grey) == 0x1);

		// m_background_count has the lower bit.
		return bool(m_background_count.load(order) & colour_mask);
	}


	typedef std::uint32_t path_length_type;


	class huddinge_environment_statistics
	{
	private:
		lb::copyable_atomic <path_length_type> m_max_path_length{};

	public:
		path_length_type path_length() const { return m_max_path_length.load(std::memory_order_relaxed); }
		path_length_type update_path_length(path_length_type length) { return update_maximum(m_max_path_length, length); }
	};


	// We assume that inserting near the (right) end of a B-tree is faster than into an arbitrary
	// position. Hence we sort by length as we process the k-mers in that order.
	struct cmp_gapmer_length_value
	{
		constexpr auto make_tuple_representation(gapmer_type gg) const { return std::make_tuple(gg.length(), gg); }

		constexpr bool operator()(gapmer_type const lhs, gapmer_type const rhs) const { return make_tuple_representation(lhs) < make_tuple_representation(rhs); }
	};


	typedef absl::btree_set <gapmer_type, cmp_gapmer_length_value> gapmer_set;
	typedef absl::btree_map <gapmer_type, count_pair, cmp_gapmer_length_value> gapmer_count_map;
	typedef absl::btree_map <gapmer_type, huddinge_environment_statistics, cmp_gapmer_length_value> gapmer_huddinge_environment_statistic_map;
	typedef absl::btree_multimap <gapmer_type, gapmer_type, cmp_gapmer_length_value> edge_map;


	struct large_count_pair
	{
	private:
		std::int64_t m_signal_count{};
		std::int64_t m_background_count{};

	public:
		void add(count_pair cc);
		std::int64_t sum() const { return m_signal_count + m_background_count; }
	};


	void large_count_pair::add(count_pair cc)
	{
		m_signal_count += cc.signal_count();
		m_background_count += cc.background_count();
	}


	large_count_pair total_count(gapmer_count_map const &mm)
	{
		large_count_pair retval{};
		for (auto const &kv : mm)
			retval.add(kv.second);

		return retval;
	}


	class count_by_length_map
	{
	private:
		typedef std::vector <large_count_pair> count_vector;

	public:
		typedef count_vector::reference reference;
		typedef count_vector::const_reference const_reference;

	private:
		count_vector m_counts;
		std::size_t m_min_k{};

	public:
		count_by_length_map(std::size_t min_k, std::size_t max_k):
			m_counts(max_k - min_k + 1),
			m_min_k{min_k}
		{
		}

		reference operator[](std::size_t ii) { return m_counts[ii - m_min_k]; }
		const_reference operator[](std::size_t ii) const { return m_counts[ii - m_min_k]; }
	};


	struct counting_context
	{
		std::size_t partition_index{};
		std::size_t read_index{};
		std::size_t lhs_position{};
	};


	sf::huddinge_distance_return_value
	huddinge_distance(
		gapmer_type gg,
		std::vector <std::uint64_t> const &expected,
		std::uint16_t expected_length
	)
	{
		std::array <std::uint64_t, 3> tested{};
		auto const tested_length{gg.length() + gg.gap_length()};
		libbio_assert_lte(tested_length, 16 * tested.size());
		gg.write_4bit_coded_to_buffer(std::span <std::uint64_t>{tested});

		return sf::huddinge_distance_4bit(
			std::span <uint64_t const>{expected},
			std::span <uint64_t const>{tested},
			expected_length,
			tested_length
		);
	}


	double huddinge_similarity(
		gapmer_type gg,
		std::vector <std::uint64_t> const &expected,
		std::uint16_t expected_length
	)
	{
		auto const res{huddinge_distance(gg, expected, expected_length)};
		return double(res.distance) / res.max_defined_characters;
	}


	// Read the k-mers and handle all in all threads (instead of partitioning).
	template <typename t_cb>
	void read_gapmers(
		sf::packed_read_vector const &reads,
		std::uint8_t const kk,
		std::size_t partition_count,
		t_cb &&cb
	)
	{
		auto const read_count{reads.size()};
		//auto const partition_size{ceil_div(read_count, partition_count)};

#pragma omp parallel for
		for (std::size_t partition_idx = 0; partition_idx < partition_count; ++partition_idx)
		{
#if 0
			auto ii{partition_idx * partition_size};
			auto const limit{std::min(read_count, (1 + partition_idx) * partition_size)};
#else
			std::size_t ii{};
			std::size_t const limit{read_count};
#endif
			for (; ii < limit; ++ii)
			{
				auto const &read{reads[ii]};
				if (read.length < kk) continue;

				counting_context ctx{partition_idx, ii};
				gapmer_type gg(read.packed_characters.front() >> (64U - 2U * kk), kk);
				cb(gg, ctx);

				read.iterate_characters(kk, [&](std::uint8_t const cc) {
					gg = gg.next_(cc);
					++ctx.lhs_position;
					cb(gg, ctx);
				});

				std::uint8_t gap_s(middle_gap_only ? kk / 2 : 1);
				auto const gap_lim{libbio::min_ct(read.length, middle_gap_only ? (kk + 3U) / 2U : kk)};

				for (; gap_s < gap_lim; ++gap_s)
				{
					for (std::uint8_t gap_l{1}; gap_l <= max_gap && kk + gap_l <= read.length; ++gap_l)
					{
						gg = gapmer_type(read.packed_characters, kk, gap_s, gap_l);
						ctx.lhs_position = 0;
						cb(gg, ctx);
						read.iterate_character_pairs(gap_s, kk + gap_l, [&](std::uint8_t const lhsc, std::uint8_t const rhsc) {
							gg = gg.next_(lhsc, rhsc);
							++ctx.lhs_position;
							cb(gg, ctx);
						});
					}
				}
			}
		}
	}


	struct seed_statistics
	{
		sf::math::p_value_type ac_test_result{};
		std::int64_t count_difference{};
		double normalised_count_difference{};
		double fold_change{};
		double similarity{};
		path_length_type path_length{};
	};


	struct current_k_tag { constexpr static bool considers_previous_k() { return false; }};
	struct prev_k_tag { constexpr static bool considers_previous_k() { return true; }};


	template <typename t_value>
	struct pearson_correlation_coefficient_accumulator
	{
		typedef t_value value_type;

		lb::accumulators::pearson_correlation_coefficient <t_value> pc{};

		virtual void operator()(seed_statistics const &statistics, first_iteration_tag) = 0;
		virtual void operator()(seed_statistics const &statistics, subsequent_iteration_tag) = 0;
		bool empty() const { return pc.empty(); }
	};


	struct count_diff_accumulator final : public pearson_correlation_coefficient_accumulator <double>
	{
		constexpr static std::string_view name() { return "Count difference v. similarity"; }

		virtual void operator()(seed_statistics const &statistics, first_iteration_tag)
		{
			pc.init(statistics.count_difference, statistics.similarity);
		}

		virtual void operator()(seed_statistics const &statistics, subsequent_iteration_tag)
		{
			pc.update(statistics.count_difference, statistics.similarity);
		}
	};


	struct normalised_count_diff_accumulator final : public pearson_correlation_coefficient_accumulator <double>
	{
		constexpr static std::string_view name() { return "Normalised count difference v. similarity"; }

		virtual void operator()(seed_statistics const &statistics, first_iteration_tag)
		{
			pc.init(statistics.normalised_count_difference, statistics.similarity);
		}

		virtual void operator()(seed_statistics const &statistics, subsequent_iteration_tag)
		{
			pc.update(statistics.normalised_count_difference, statistics.similarity);
		}
	};


	struct fold_change_accumulator final : public pearson_correlation_coefficient_accumulator <double>
	{
		constexpr static std::string_view name() { return "Fold change v. similarity"; }

		virtual void operator()(seed_statistics const &statistics, first_iteration_tag)
		{
			pc.init(statistics.fold_change, statistics.similarity);
		}

		virtual void operator()(seed_statistics const &statistics, subsequent_iteration_tag)
		{
			pc.update(statistics.fold_change, statistics.similarity);
		}
	};


	struct ac_test_accumulator final : public pearson_correlation_coefficient_accumulator <sf::math::p_value_type>
	{
		constexpr static std::string_view name() { return "Audic-Claverie test v. similarity"; }

		virtual void operator()(seed_statistics const &statistics, first_iteration_tag)
		{
			pc.init(statistics.ac_test_result, statistics.similarity);
		}

		virtual void operator()(seed_statistics const &statistics, subsequent_iteration_tag)
		{
			pc.update(statistics.ac_test_result, statistics.similarity);
		}
	};


	struct path_length_accumulator final : public pearson_correlation_coefficient_accumulator <sf::math::p_value_type>
	{
		constexpr static std::string_view name() { return "Path length v. similarity"; }

		virtual void operator()(seed_statistics const &statistics, first_iteration_tag)
		{
			pc.init(statistics.path_length, statistics.similarity);
		}

		virtual void operator()(seed_statistics const &statistics, subsequent_iteration_tag)
		{
			pc.update(statistics.path_length, statistics.similarity);
		}
	};


	typedef std::tuple <
		count_diff_accumulator,
		normalised_count_diff_accumulator,
		fold_change_accumulator,
		ac_test_accumulator
	> accumulator_tuple_;

	typedef lb::tuples::append_t <
		accumulator_tuple_,
		path_length_accumulator
	> semi_local_accumulator_tuple_;


	template <typename t_accumulators>
	struct accumulator_tuple_tpl
	{
		t_accumulators accumulators{};
		t_accumulators &tie() { return accumulators; }
		t_accumulators const &tie() const { return accumulators; }

		bool empty() const { return std::get <0>(accumulators).empty(); }

		template <typename t_iteration_tag>
		void update(seed_statistics const &statistics, t_iteration_tag iteration_tag);

		void init_or_update(seed_statistics const &statistics);

		void clear() { accumulators = t_accumulators{}; }
	};


	template <typename t_accumulators>
	struct named_accumulator_tuple_tpl : public accumulator_tuple_tpl <t_accumulators>
	{
		std::string name;

		explicit named_accumulator_tuple_tpl(std::string_view name_):
			name{name_}
		{
		}

		void output(std::ostream &os, std::uint8_t kk) const;
	};


	typedef accumulator_tuple_tpl <accumulator_tuple_> accumulator_tuple;


	template <typename t_accumulators>
	template <typename t_iteration_tag>
	void accumulator_tuple_tpl <t_accumulators>::update(seed_statistics const &statistics, t_iteration_tag iteration_tag)
	{
		// Update the accumulators.
		std::apply([&](auto && ... acc) -> void {
			(acc(statistics, t_iteration_tag{}), ...);
		}, tie());
	}


	template <typename t_accumulators>
	void accumulator_tuple_tpl <t_accumulators>::init_or_update(seed_statistics const &statistics)
	{
		if (empty())
		{
			std::apply([&](auto && ... acc) -> void {
				(acc(statistics, first_iteration_tag{}), ...);
			}, tie());
		}
		else
		{
			std::apply([&](auto && ... acc) -> void {
				(acc(statistics, subsequent_iteration_tag{}), ...);
			}, tie());
		}
	}


	template <typename t_accumulators>
	void named_accumulator_tuple_tpl <t_accumulators>::output(std::ostream &os, std::uint8_t kk) const
	{
		std::apply([&](auto && ... acc) -> void {
			((std::cout << acc.name() << '\t' << name << '\t' << +kk << '\t' << acc.pc.value() << '\n'), ...);
		}, this->tie());
	}


	struct accumulator_group
	{
		typedef named_accumulator_tuple_tpl <accumulator_tuple_> accumulator_tuple_type;
		typedef named_accumulator_tuple_tpl <semi_local_accumulator_tuple_> semi_local_accumulator_tuple_type;

		accumulator_tuple_type all_nodes;
		//accumulator_tuple_type all_nodes_current_k; // XXX This yields nonsensical results since we do not correct the effect of the similarities depending on k.
		accumulator_tuple_type local_maximum;
		semi_local_accumulator_tuple_type semi_local_maximum;

		accumulator_group():
			all_nodes{"All"},
			//all_nodes_current_k{"All (current k)"},
			local_maximum{"Local maximum"},
			semi_local_maximum{"Semi-local maximum"}
		{
		}

		void reset();
		void output(std::ostream &os, std::uint8_t kk) const;

	private:
		decltype(auto) tie() const { return std::tie(all_nodes, /* all_nodes_current_k, */ local_maximum, semi_local_maximum); }
	};


	void accumulator_group::reset()
	{
		//all_nodes_current_k.clear();

		// Local and semi-local maximum statistics are re-calculated on each iteration.
		local_maximum.clear();
		semi_local_maximum.clear();
	}


	void accumulator_group::output(std::ostream &os, std::uint8_t kk) const
	{
		std::apply([&](auto && ... acc) -> void {
			(acc.output(std::cout, kk), ...);
		}, tie());
	}


	void count_gapmers(
		sf::packed_read_vector const &signal_reads,
		sf::packed_read_vector const &background_reads,
		std::uint8_t const kk,
		gapmer_count_map &counts
	)
	{
		typedef std::vector <gapmer_count_map> buffer_vector;
		buffer_vector buffers(16);

		auto const do_process([kk, &buffers](auto const &reads, auto &&increment){
			read_gapmers(reads, kk, buffers.size(), [&](gapmer_type const gg, counting_context const &ctx){
				// Partition by k-mer suffix.
				// FIXME: Partitioning by prefix would probably be a better idea but due to gapmer’s layout it is also more difficult.
				auto const key{gg.value() & UINT64_C(0xf)};
				libbio_assert_lt(key, buffers.size());
				if (key == ctx.partition_index)
				{
					auto &buffer{buffers[key]};
					increment(gg, buffer);
				}
			});
		});

		do_process(signal_reads, [](gapmer_type gg, gapmer_count_map &mm){ mm[gg].increment_signal(); });
		do_process(background_reads, [](gapmer_type gg, gapmer_count_map &mm){ mm[gg].increment_background(); });

		while (!buffers.empty())
		{
			counts.merge(buffers.back());
			buffers.pop_back();
		}
	}


	// Error reporter for FASTA/Q reader adapter.
	struct reader_adapter_delegate : public sf::reader_adapter_delegate
	{
		bool should_report_errors_for_path(sf::reader_adapter &, std::string_view path) override
		{
			return true;
		}

		void found_first_read_with_unexpected_character(
			sf::reader_adapter&, std::string_view path,
			std::uint64_t lineno
		) override
		{
			std::cerr << "WARNING: Skipping reads with unexpected characters in " << path << "; first one on line " << lineno << ".\n";
		}

		void found_total_reads_with_unexpected_characters(
			sf::reader_adapter&,
			std::string_view path,
			std::uint64_t count
		) override
		{
			std::cerr << "WARNING: Skipped " << count << " reads in " << path << ".\n";
		}
	};


	struct configuration
	{
		std::string signal_path;
		std::string background_path;
		std::string tested_kmer;
		std::string graphviz_output_prefix;
		std::uint16_t max_k{gapmer_type::max_k}; // args.hxx outputs this in a nice way if not std::uint8_t.
	};


	configuration process_command_line_options(int argc, char **argv)
	{
		configuration retval;
		args::ArgumentParser parser(
			std::format(
				"Calculate the correlations of Huddinge similarity to the expected gapmer and "
				"various properties. Intermediate values will be written to stdout during the "
				"calculation. Max " SF_FORMAT_ITALIC "k" SF_FORMAT_RESET " = {}, max gap size "
				"= {}, middle gap only = {}.\n\n"
				"We say that a node is a " SF_FORMAT_ITALIC "local maximum" SF_FORMAT_RESET
				" of its Huddinge-1 neighbourhood (up to the current " SF_FORMAT_ITALIC "k"
				SF_FORMAT_RESET ") if the adjacent nodes have lower fold change.\n\n"
				"We say that a node " SF_FORMAT_ITALIC "v(i(l))" SF_FORMAT_RESET " is a "
				SF_FORMAT_ITALIC "semi-local maximum" SF_FORMAT_RESET " of its Huddinge-1 "
				"neighbourhood (up to the current " SF_FORMAT_ITALIC "k" SF_FORMAT_RESET ") if "
				"there exists a path " SF_FORMAT_ITALIC "v(i(1)), … v(i(l))" SF_FORMAT_RESET
			 	" such that the corresponding sequence of values (" SF_FORMAT_ITALIC "f"
				SF_FORMAT_RESET ", " SF_FORMAT_ITALIC "k’" SF_FORMAT_RESET ") where "
				SF_FORMAT_ITALIC "f" SF_FORMAT_RESET " is the fold change and " SF_FORMAT_ITALIC
				"k’" SF_FORMAT_RESET " is the corresponding value of " SF_FORMAT_ITALIC "k"
				SF_FORMAT_RESET " is strictly increasing (compared left to right) and that |"
				SF_FORMAT_ITALIC "k’ⱼ" SF_FORMAT_RESET "₋₁ - " SF_FORMAT_ITALIC "k’ⱼ"
				SF_FORMAT_RESET "| ≤ 1 for all " SF_FORMAT_ITALIC "j" SF_FORMAT_RESET ".",
				gapmer_type::max_k,
				gapmer_type::max_gap,
				gapmer_type::middle_gap_only
			)
		);
		parser.SetArgumentSeparations(true, true, true, true);

		{
			args::HelpParams help_params{};
			help_params.addDefault = true;
			help_params.shortPrefix = "-";
			help_params.longPrefix = "--";
			help_params.shortSeparator = " ";
			help_params.longSeparator = " ";
			help_params.defaultString = " Default: ";
			parser.helpParams = help_params;
		}

		args::HelpFlag help_(parser, "help", "Display this help.", {'h', "help"});
		args::CompletionFlag completion_(parser, {"complete"});

		sf::args::value_flag background_path_(
			parser,
			"path",
			"Background FASTA file.",
			{'b', "background"},
			retval.background_path,
			args::Options::Required
		);

		sf::args::value_flag signal_path_(
			parser,
			"path",
			"Signal FASTA file.",
			{'s', "signal"},
			retval.signal_path,
			args::Options::Required
		);

		sf::args::value_flag tested_kmer_(
			parser,
			"k-mer",
			"k-mer to test against.",
			{'t', "tested-kmer"},
			retval.tested_kmer,
			args::Options::Required
		);

		args::ValueFlag <std::uint16_t>max_k_(
			parser,
			"max-k",
			"Maximum k to test. Default is the number of defined characters in expected plus two.",
			{'m', "max-k"},
			args::Options::None
		);
		max_k_.HelpDefault("");

		sf::args::value_flag graphviz_output_prefix_(
			parser,
			"graphviz-output-prefix",
			"Compressed Graphviz output prefix",
			{"graphviz-output-prefix"},
			retval.graphviz_output_prefix,
			args::Options::None
		);

		// Parse and check.
		try
		{
			parser.ParseCLI(argc, argv);
		}
		catch (args::Help const&)
		{
			std::cerr << parser;
			std::exit(0);
		}
		catch (args::Completion const&)
		{
			std::cerr << parser;
			std::exit(0);
		}
		catch (sf::args::output_version const&)
		{
			std::cout << "test_counting " << sf::version << '\n';
			std::exit(0);
		}
		catch (std::runtime_error const& err)
		{
			std::cerr << err.what() << '\n';
			std::exit(1);
		}

		if (max_k_)
			retval.max_k = args::get(max_k_);
		else
		{
			retval.max_k = 2 + std::count_if(
				retval.tested_kmer.begin(),
				retval.tested_kmer.end(),
				[](auto const cc){
					switch (cc)
					{
						case '-':
						case '.':
						case 'n':
							return false;
						default:
							return true;
					}
				}
			);
		}

		return retval;
	}


	// Graph classes for BGL to pass to topological sort.
	struct edge_tuple
	{
		edge_map out_edges_k_to_k;
		edge_map out_edges_k_to_k1;
		edge_map out_edges_k1_to_k;
		edge_map out_edges_k1_to_k1;

		auto tie() { return std::tie(out_edges_k_to_k, out_edges_k_to_k1, out_edges_k1_to_k, out_edges_k1_to_k1); }
		void clear() { std::apply([&](auto && ... edges){ (edges.clear(), ...); }, tie()); }
		void merge_unique(edge_tuple &src);
		void next_k();
	};


	void edge_tuple::merge_unique(edge_tuple &src)
	{
		lb::pairwise_apply(src.tie(), tie(), [&](edge_map &src, edge_map &dst){
			auto src_it{src.begin()};
			while (src_it != src.end()) // Do not cache end() since erase() can invalidate it.
			{
				// Check if the source value already exists in dst.
				auto dst_range{dst.equal_range(src_it->first)};
				auto dst_it{dst_range.first};
				auto const dst_end{dst_range.second};
				for (; dst_it != dst_end; ++dst_it)
				{
					if (*src_it == *dst_it)
					{
						// Found a match; erase from src.
						src_it = src.erase(src_it);
						goto continue_loop;
					}
				}

				// No match found.
				++src_it;

			continue_loop:
				;
			}

			// Duplicates have been handled; merge.
			dst.merge(src);
		});
	}


	void edge_tuple::next_k()
	{
		using std::swap;
		swap(out_edges_k_to_k, out_edges_k1_to_k1);

		out_edges_k_to_k.clear();
		out_edges_k_to_k1.clear();
		out_edges_k1_to_k.clear();
	}


	class h1_graph
	{
	private:
		typedef std::array <gapmer_count_map *, 2> count_map_array;
		typedef std::array <edge_map const *, 4> edge_map_array;

		count_map_array m_counts;
		edge_map_array m_edges;
		std::uint16_t m_kk{};

	public:
		typedef gapmer_type vertex_type;
		typedef std::pair <gapmer_type, gapmer_type> edge_type;

		typedef lb::indirect_join_iterator <
			count_map_array::const_iterator
		> joined_count_iterator;

		typedef lb::joined_range_iterator <
			edge_map::const_iterator,
			edge_map::const_iterator,
			std::tuple_size_v <edge_map_array> // std::array is tuple-like.
		> joined_edge_iterator;

		typedef boost::transform_iterator <
			access_map_key,
			joined_count_iterator
		> node_iterator;

		typedef std::pair <node_iterator, node_iterator> node_iterator_pair;
		typedef std::pair <joined_edge_iterator, joined_edge_iterator> joined_edge_iterator_pair;

	private:
		constexpr static bool is_current_k_source_edge_array_index(std::size_t idx) { return idx < 2; }

	public:
		h1_graph(
			gapmer_count_map &counts,
			gapmer_count_map &prev_counts,
			edge_tuple const &edges
		):
			m_counts{&counts, &prev_counts},
			m_edges{ // is_current_k_source_edge_array_index() depends on order.
				&edges.out_edges_k_to_k,
				&edges.out_edges_k_to_k1,
				&edges.out_edges_k1_to_k,
				&edges.out_edges_k1_to_k1
			}
		{
		}

		void set_k(std::uint16_t kk) { m_kk = kk; }

		std::size_t node_count() const;
		std::size_t edge_count() const;
		std::size_t out_degree(gapmer_type gg) const;
		node_iterator_pair node_range() const;
		joined_edge_iterator_pair edge_range() const;
		joined_edge_iterator_pair out_edges(gapmer_type gg) const;

		count_pair &find_count(gapmer_type gg);
	};


	struct h1_graph_colour_map
	{
		typedef gapmer_type key_type;
		typedef colour_type value_type;
		typedef value_type reference;
		typedef boost::read_write_property_map_tag category;

		h1_graph *graph{};

		explicit h1_graph_colour_map(h1_graph &graph_): graph{&graph_} {}

		inline value_type get(key_type gg) const;
		inline void put(key_type gg, value_type cc);
	};


	struct h1_graph_vertex_id_map
	{
		typedef gapmer_type key_type;
		typedef std::uint64_t value_type;
		typedef value_type reference;
		typedef boost::readable_property_map_tag category;

		inline value_type get(key_type gg) const { return std::uint64_t(gg); };
	};


	class h1_graph_property_writer
	{
	private:
		h1_graph *m_graph{};
		std::uint64_t m_signal_size{};
		std::uint64_t m_background_size{};

	public:
		h1_graph_property_writer(h1_graph &graph, std::uint64_t signal_size, std::uint64_t background_size):
			m_graph{&graph},
			m_signal_size{signal_size},
			m_background_size{background_size}
		{
		}

		void operator()(std::ostream &os) const;
		void operator()(std::ostream &os, h1_graph::vertex_type gapmer) const;
		void operator()(std::ostream &os, h1_graph::edge_type gapmer_pair) const;
	};


	count_pair &h1_graph::find_count(gapmer_type gg)
	{
		auto const kk{gg.length()};
		libbio_assert(kk == m_kk or kk + 1 == m_kk);

		auto const do_find{[&](auto &counts) -> count_pair & {
			auto const it{counts.find(gg)};
			libbio_assert_neq(it, counts.end());
			return it->second;
		}};

		if (kk == m_kk)
			return do_find(*std::get <0>(m_counts)); // current
		else
			return do_find(*std::get <1>(m_counts)); // prev
	}


	std::size_t h1_graph::node_count() const
	{
		return std::accumulate(m_counts.begin(), m_counts.end(), std::size_t{}, [](std::size_t acc, auto const &mm){
			return acc + mm->size();
		});
	}


	std::size_t h1_graph::edge_count() const
	{
		return std::accumulate(m_edges.begin(), m_edges.end(), std::size_t{}, [](std::size_t acc, auto const &mm){
			return acc + mm->size();
		});
	}


	std::size_t h1_graph::out_degree(gapmer_type gg) const
	{
		auto const kk_{gg.length()};
		libbio_assert(m_kk == kk_ or m_kk == kk_ + 1);
		bool const should_use_current{m_kk == kk_};

		std::size_t idx{};
		return std::accumulate(m_edges.begin(), m_edges.end(), std::size_t{}, [&](std::size_t acc, auto const &mm){
			++idx;
			if (should_use_current == is_current_k_source_edge_array_index(idx - 1))
				return acc + mm->count(gg);
			return acc;
		});
	}


	auto h1_graph::node_range() const -> node_iterator_pair
	{
		auto iterators{lb::make_indirect_join_iterator_pair(m_counts)};
		return {
			node_iterator{iterators.first, access_map_key{}},
			node_iterator{iterators.second, access_map_key{}}
		};
	}


	auto h1_graph::edge_range() const -> joined_edge_iterator_pair
	{
		auto begins{std::apply([](auto && ... mm){ return std::array{mm->cbegin()...}; }, m_edges)};
		auto ends{std::apply([](auto && ... mm){ return std::array{mm->cend()...}; }, m_edges)};
		joined_edge_iterator it{std::move(begins), std::move(ends)};
		auto end{it.make_sentinel()};
		return {std::move(it), std::move(end)};
	}


	auto h1_graph::out_edges(gapmer_type gg) const -> joined_edge_iterator_pair
	{
		auto const kk_{gg.length()};
		libbio_assert(m_kk == kk_ or m_kk == kk_ + 1);
		bool const should_use_current{m_kk == kk_};
		auto begins{lb::map_to_array(m_edges, [should_use_current, gg](auto const idx, edge_map_array::value_type mm){
			if (should_use_current == is_current_k_source_edge_array_index(idx))
				return mm->lower_bound(gg);
			else
				return mm->end();
		})};
		auto ends{lb::map_to_array(m_edges, [should_use_current, gg](auto const idx, edge_map_array::value_type mm){
			if (should_use_current == is_current_k_source_edge_array_index(idx))
				return mm->upper_bound(gg);
			else
				return mm->end();
		})};
		joined_edge_iterator it{std::move(begins), std::move(ends)};
		auto end{it.make_sentinel()};
		return {std::move(it), std::move(end)};
	}


	auto h1_graph_colour_map::get(key_type gg) const -> value_type
	{
		return graph->find_count(gg).colour();
	}


	void h1_graph_colour_map::put(key_type gg, value_type cc)
	{
		graph->find_count(gg).set_colour(cc);
	}


	void h1_graph_property_writer::operator()(std::ostream &os) const
	{
	}


	void h1_graph_property_writer::operator()(std::ostream &os, h1_graph::vertex_type gg) const
	{
		auto const count{m_graph->find_count(gg)};
		auto const fold_change{count.fold_change(m_signal_size, m_background_size)};

		os << std::format(
			"[label = \"{:x}\\nSC: {}\\nBC: {}\\nFC: {}\"]",
			std::uint64_t(gg),
			count.signal_count(),
			count.background_count(),
			fold_change
		);
	}


	void h1_graph_property_writer::operator()(std::ostream &os, h1_graph::edge_type gapmer_pair) const
	{
	}


	// Boost concepts; functions found via ADL.
	// VertexListGraph
	h1_graph::node_iterator_pair vertices(h1_graph const &graph) { return graph.node_range(); }
	std::size_t num_vertices(h1_graph const &graph) { return graph.node_count(); }

	// EdgeListGraph
	h1_graph::joined_edge_iterator_pair edges(h1_graph const &graph) { return graph.edge_range(); }
	std::size_t num_edges(h1_graph const &graph) { return graph.edge_count(); }

	// IncidenceGraph
	gapmer_type source(edge_map::value_type edge, h1_graph const &) { return edge.first; }
	gapmer_type target(edge_map::value_type edge, h1_graph const &) { return edge.second; }
	h1_graph::joined_edge_iterator_pair out_edges(gapmer_type gg, h1_graph const &graph) { return graph.out_edges(gg); }
	std::size_t out_degree(gapmer_type gg, h1_graph const &graph) { return graph.out_degree(gg); }

	// Properties
	inline h1_graph_vertex_id_map::value_type get(
		h1_graph_vertex_id_map const &map,
		h1_graph_vertex_id_map::key_type key
	)
	{
		return map.get(key);
	}

	inline h1_graph_colour_map::value_type get(
		h1_graph_colour_map const &map,
		h1_graph_colour_map::key_type key
	)
	{
		return map.get(key);
	}

	inline void put(
		h1_graph_colour_map &map,
		h1_graph_colour_map::key_type key,
		colour_type value
	)
	{
		map.put(key, value);
	}


	// For debugging.
	struct cycle_detector : public boost::dfs_visitor <>
	{
		bool has_cycle{};

		template <typename t_edge, typename t_graph>
		void back_edge(t_edge const &edge, t_graph const &graph)
		{
			has_cycle = true;
			std::cerr << std::format("ERROR: Cycle detected: {:x} → {:x}\n", std::uint64_t(edge.first), std::uint64_t(edge.second));
		}
	};


	// For parallelising.
	struct partition_meta
	{
		accumulator_tuple accumulators;
		gapmer_vector found_seeds_current_k;
		gapmer_vector found_seeds_prev_k;
		gapmer_vector discarded_seeds;
		edge_tuple edges;

		void reset()
		{
			accumulators = accumulator_tuple{};
			found_seeds_prev_k.clear();
			found_seeds_current_k.clear();
			discarded_seeds.clear();
			edges.clear();
		}
	};
}


namespace boost {

	// Graph traits, need to be in namespace boost.
	template <>
	struct graph_traits <h1_graph>
	{
		// For topological sorting, we need to model Graph, IncidenceGraph and VertexListGraph.
		// We also model EdgeListGraph to get Graphviz output.
		// IncidenceGraph requires the source() function and hence we represent edges with a
		// pair of nodes.
		struct traversal_category : public vertex_list_graph_tag, public edge_list_graph_tag, public incidence_graph_tag {};
		typedef h1_graph::vertex_type vertex_descriptor;
		typedef h1_graph::edge_type edge_descriptor;

		typedef directed_tag directed_category;
		typedef disallow_parallel_edge_tag edge_parallel_category;
		typedef std::size_t vertices_size_type;
		typedef std::size_t edges_size_type;

		typedef h1_graph::node_iterator vertex_iterator;
		typedef h1_graph::joined_edge_iterator edge_iterator;
		typedef h1_graph::joined_edge_iterator out_edge_iterator;
		typedef std::size_t degree_size_type;

		constexpr static vertex_descriptor null_vertex() { return vertex_descriptor::make_marked(); }
	};

	// Colour traits.
	template <>
	struct color_traits <colour_type>
	{
		static colour_type white() { return colour_type::white; }
		static colour_type gray() { return colour_type::grey; }
		static colour_type black() { return colour_type::black; }
	};
}


int main(int argc, char **argv)
{
	std::ios_base::sync_with_stdio(false);

	auto const configuration{process_command_line_options(argc, argv)};

	sf::packed_read_vector signal_reads;
	sf::packed_read_vector background_reads;

	// Process the input reads.
	lb::log_time(std::cerr) << "Loading input reads…\n";
	{
		reader_adapter_delegate delegate;
		sf::reader_adapter_type reader{delegate};

		auto const process_path{
			[&](std::string const& path, sf::packed_read_vector& dst) {
				sf::reader_adapter_guard guard{reader};
				reader.read_from_path(path);
				while (reader.retrieve_next_read())
					dst.emplace_back(reader.read_buffer(), reader.read_length());
			}
		};

		process_path(configuration.signal_path, signal_reads);
		process_path(configuration.background_path, background_reads);
	}

	auto const read_length_sum{[](sf::packed_read_vector const& reads) -> std::uint64_t {
		return std::accumulate(
			reads.begin(),
			reads.end(),
			std::uint64_t{},
			[](auto const acc, sf::packed_read const& rr) {
				return acc + rr.length;
			}
		);
	}};

	auto const signal_size{read_length_sum(signal_reads)};
	auto const background_size{read_length_sum(background_reads)};
	auto const signal_to_total_length_ratio{double(signal_size) / (signal_size + background_size)};

	if (count_pair::max_count() < signal_size or count_pair::max_count() < background_size)
	{
		std::cerr << "WARNING: Either uncompressed signal or background length exceeds maximum k-mer count ("
			<< signal_size
			<< ", "
			<< background_size
			<< ", "
			<< count_pair::max_count()
			<< "); overflows are possible.\n";
	}

	auto const expected{[&]{
		packed_word_vector expected{};
		sf::pack_characters_ <sf::dna_alphabet::dna16>(configuration.tested_kmer, expected);
		return expected;
	}()};
	auto const expected_length{configuration.tested_kmer.size()};

	constexpr std::uint8_t min_k{5};
	auto const max_k{lb::min_ct(gapmer_type::max_k, configuration.max_k)};
	lb::log_time(std::cerr) << "Using " << +max_k << " for maximum k.\n";

	gapmer_count_map counts;
	gapmer_count_map prev_counts;
	accumulator_group accumulators;
	gapmer_set local_maxima;
	gapmer_huddinge_environment_statistic_map semi_local_maximum_statistics;
	count_by_length_map total_counts{min_k, max_k};
	edge_tuple edges;
	h1_graph graph{counts, prev_counts, edges};
	gapmer_vector sorted_gapmers;
	gapmer_vector found_semi_local_maxima;
	std::mutex found_semi_local_maxima_mutex; // FIXME: align to cache line size?
	std::atomic_uint64_t nodes_without_neighbours{};

	auto const do_count{[&](std::uint8_t const kk){
		count_gapmers(signal_reads, background_reads, kk, counts);
		total_counts[kk] = total_count(counts);
		lb::log_time(std::cerr) << "  Got " << counts.size() << " distinct k-mers.\n";
	}};

	auto const add_edge{[](edge_map &map, gapmer_type src, gapmer_type dst){
		auto pp{map.equal_range(src)};
		auto it{pp.first};
		auto const end{pp.second};

		// Check if we already have the edge.
		for (; it != end; ++it)
		{
			if (it->second == dst)
				return;
		}

		map.emplace_hint(end, src, dst);
	}};

	auto const calculate_seed_statistics{
		[&](
			gapmer_type const gg,
			count_pair const counts,
			std::uint8_t const kk
		) -> seed_statistics
		{
			seed_statistics statistics{};

			// Audic-Claverie
			{
				auto const counts_{counts.with_pseudocounts()};
				auto const res{
					sf::math::beta_incomplete(
						counts_.signal_count(),
						counts_.background_count(),
						signal_to_total_length_ratio
					)
				};

				if (res)
					statistics.ac_test_result = res.value;
				else
				{
					std::cerr << "Got error " << res.error << " while doing A-C test.\n";
					statistics.ac_test_result = std::numeric_limits <sf::math::p_value_type>::max();
				}
			}

			// Count diffs
			statistics.count_difference = counts.diff();
			statistics.normalised_count_difference = double(counts.diff()) / total_counts[kk].sum();

			// Fold change
			statistics.fold_change = counts.fold_change(signal_size, background_size);

			// Similarity to the expected k-mer
			statistics.similarity = huddinge_similarity(gg, expected, expected_length);

			return statistics;
		}
	};

	// Count the k-mers for each k and process.
	{
		auto const hardware_concurrency{std::thread::hardware_concurrency()};
		auto const partition_count{2U * hardware_concurrency};

		std::vector <partition_meta> partition_specifics(partition_count);

		lb::log_time(std::cerr) << "Counting k-mers…\n";
		std::cout << "Accumulator\tNodes considered\tk\tCorrelation\n";
		lb::log_time(std::cerr) << " Handling k = " << +min_k << '\n';
		do_count(min_k);
		for (std::uint8_t kk{min_k + 1U}; kk <= max_k; ++kk)
		{
			lb::log_time(std::cerr) << " Handling k = " << +kk << '\n';

			// Reset the state.
			{
				using std::swap;
				swap(counts, prev_counts);

				counts.clear();
				sorted_gapmers.clear();

				accumulators.reset();

				edges.next_k();
				graph.set_k(kk);
			}

			// Reset the partition-specific data.
			for (auto &meta : partition_specifics)
				meta.reset();

			do_count(kk);

			auto const distinct_gapmer_count{counts.size()};
			auto const partition_size{ceil_div(distinct_gapmer_count, partition_count)};

			// Process each partition.
			lb::log_time(std::cerr) << "  Processing k-mers…\n";
#pragma omp parallel for
			for (std::size_t partition_idx = 0; partition_idx < partition_count; ++partition_idx)
			{
				auto const start_idx{partition_idx * partition_size};
				auto const end_idx{std::min(distinct_gapmer_count, start_idx + partition_size)};

				// Get an iterator pair for the current partition.
				// Using operator+= should be reasonably fast in Abseil.
				auto it{counts.begin()};
				it += start_idx;
				auto const end{[&]{
					auto retval{counts.begin()};
					retval += end_idx;
					return retval;
				}()};

				libbio_assert_lt(partition_idx, partition_specifics.size());
				auto &spec{partition_specifics[partition_idx]};

				for_each_(it, end, [&](auto it, auto const iteration_tag){
					auto const current_gapmer{it->first};
					auto const current_counts{it->second};
					auto const statistics{calculate_seed_statistics(current_gapmer, current_counts, kk)};
					spec.accumulators.update(statistics, iteration_tag);

					// Huddinge-1 environment
					// Since we now have the counts, we can calculate the enrichment values.
					bool has_neighbours{};
					auto const handle_neighbour{
						[&](
							gapmer_type const neighbour,
							edge_map &current_out_edges,
							edge_map &neighbour_out_edges,
							gapmer_vector &found_neighbour_seeds,
							auto const k_tag
						) -> void {
							auto const it{counts.find(neighbour)};
							if (counts.end() != it)
							{
								has_neighbours = true;

								// Since determining the size of the Huddinge-1 neighbourhood is difficult while saving memory,
								// we skip normalising the fold change values.
								auto const neighbour_fold_change{it->second.fold_change(signal_size, background_size)};

								// Check if the current gapmer’s fold change is greater than neighbour’s.
								if (neighbour_fold_change == statistics.fold_change)
								{
									if constexpr (k_tag.considers_previous_k())
									{
										// If we are comparing to the previous k, add the edge.
										add_edge(current_out_edges, neighbour, current_gapmer);
										if (local_maxima.contains(neighbour))
											spec.discarded_seeds.push_back(neighbour);
									}
								}
								else if (neighbour_fold_change < statistics.fold_change)
								{
									add_edge(current_out_edges, neighbour, current_gapmer);
									if constexpr (k_tag.considers_previous_k())
									{
										if (local_maxima.contains(neighbour))
											spec.discarded_seeds.push_back(neighbour);
									}
								}
								else
								{
									add_edge(neighbour_out_edges, current_gapmer, neighbour);
									found_neighbour_seeds.push_back(neighbour);
									spec.discarded_seeds.push_back(current_gapmer);
								}
							}
						}
					};

					// Edges between k and k - 1
					current_gapmer.template huddinge_neighbours <false, true, true>([&](gapmer_type neighbour){
						handle_neighbour(
							neighbour,
							spec.edges.out_edges_k_to_k1,
							spec.edges.out_edges_k1_to_k,
							spec.found_seeds_prev_k,
							prev_k_tag{}
						);
					});

					// Edges for k
					current_gapmer.template huddinge_neighbours <true, false, true>([&](gapmer_type neighbour){
						handle_neighbour(
							neighbour,
							spec.edges.out_edges_k_to_k,
							spec.edges.out_edges_k_to_k,
							spec.found_seeds_current_k,
							current_k_tag{}
						);
					});

					if (not has_neighbours)
						nodes_without_neighbours.fetch_add(1, std::memory_order_relaxed);

					// Sort to optimise handling in the main thread.
					std::sort(spec.found_seeds_prev_k.begin(), spec.found_seeds_prev_k.end(), cmp_gapmer_length_value{});
					std::sort(spec.found_seeds_current_k.begin(), spec.found_seeds_current_k.end(), cmp_gapmer_length_value{});
					std::sort(spec.discarded_seeds.begin(), spec.discarded_seeds.end(), cmp_gapmer_length_value{});
				});
			}

			lb::log_time(std::cerr) << "  Total number of nodes without neighbours: " << nodes_without_neighbours << '\n';
			lb::log_time(std::cerr) << "  Removing false positives from found seeds…\n";
			auto const make_check_function{[](auto &seeds){
				return [it = seeds.begin(), end = seeds.end()](auto const discarded_seed) mutable {
					libbio_assert(not discarded_seed.is_marked());
					it = std::lower_bound(it, end, discarded_seed, [](gapmer_type lhs, gapmer_type const rhs){
						// rhs should always be discarded_seed but we play safe.
						return cmp_gapmer_length_value{}(lhs.with_mark_cleared(), rhs.with_mark_cleared());
					});

					while (it != end && it->with_mark_cleared() == discarded_seed)
					{
						it->mark();
						++it;
					}
				};
			}};

#pragma omp parallel for
			for (std::size_t partition_idx = 0; partition_idx < partition_count; ++partition_idx)
			{
				// Check the discarded seeds in all partitions and remove from the current one’s found seeds.
				// This is probably faster than first inserting and then erasing,
				// since modifying the B-tree is not particularly fast.

				auto &spec{partition_specifics[partition_idx]};
				auto &found_seeds{spec.found_seeds_current_k};
				auto &prev_found_seeds{spec.found_seeds_prev_k};

				for (std::size_t partition_idx_ = 0; partition_idx_ < partition_count; ++partition_idx_)
				{
					// Re-instantiate every time since we start the search from the beginning.
					auto check_fn{make_check_function(found_seeds)};
					auto prev_check_fn{make_check_function(prev_found_seeds)};

					auto const &discarded_seeds{partition_specifics[partition_idx_].discarded_seeds};
					for (auto const gg : discarded_seeds)
					{
						if (kk == gg.length())
							check_fn(gg);
						else
						{
							libbio_assert_eq(gg.length() + 1, +kk);
							prev_check_fn(gg);
						}
					}
				}
			}

			// Update using the partition data.
			lb::log_time(std::cerr) << "  Merging collected data…\n";
			for (auto &spec : partition_specifics)
			{
				// Keep the non-marked gapmers.
				auto const filter_fn{[](gapmer_type const gg){ return not gg.is_marked(); }};

				// Handle the local maxima.
				// Note that spec.found_seeds and spec.discarded_seeds need to be sorted.
				for (auto const gg : spec.discarded_seeds)
					local_maxima.erase(gg);

				for (auto const gg : spec.found_seeds_prev_k | rsv::filter(filter_fn))
					local_maxima.insert(gg);

				for (auto const gg : spec.found_seeds_current_k | rsv::filter(filter_fn))
					local_maxima.insert(gg);

				// Make sure that the candidate seeds have counts.
				// (We require this in the semi-local maximum handling algorithm.)
				for (auto const gg : spec.found_seeds_prev_k | rsv::filter(filter_fn))
					prev_counts.try_emplace(gg, count_pair{});

				for (auto const gg : spec.found_seeds_current_k | rsv::filter(filter_fn))
					counts.try_emplace(gg, count_pair{});

				// These do not deallocate memory.
				spec.found_seeds_prev_k.clear();
				spec.found_seeds_current_k.clear();
				spec.discarded_seeds.clear();

				// Handle the edges.
				edges.merge_unique(spec.edges);

				// Handle the accumulators for all nodes.
				auto const do_update{[&](auto &dst){
					lb::pairwise_apply(
						spec.accumulators.tie(),
						dst.tie(),
						[](auto &src, auto &dst){
							auto &src_{src.pc};
							auto &dst_{dst.pc};

							if (src_.empty()) return;

							// src_ not empty.
							if (dst_.empty())
							{
								dst_ = src_;
								return;
							}

							// Neither src_ nor dst_ empty.
							dst_.pairwise_update(src_);
						}
					);
				}};

				do_update(accumulators.all_nodes);
				//do_update(accumulators.all_nodes_current_k);
			}

			lb::log_time(std::cerr) << "  Found " << local_maxima.size() << " local maxima.\n";

			// The graph is now valid b.c. the counted k-mers (nodes) have
			// been updated and the edge B-trees have been merged.
			if (not configuration.graphviz_output_prefix.empty())
			{
				lb::log_time(std::cerr) << "  Outputting for Graphviz…\n";
				std::string const path{std::format("{}.k{}.dot.gz", configuration.graphviz_output_prefix, +kk)};

				lb::file_ostream file;
				boost::iostreams::gzip_compressor compressor;
				boost::iostreams::filtering_ostream os;
				os.push(compressor);
				os.push(file);

				h1_graph_property_writer property_writer{graph, signal_size, background_size};
				h1_graph_vertex_id_map vertex_id_mapper{};

				lb::open_file_for_writing(path, file, lb::make_writing_open_mode({lb::writing_open_mode::CREATE}));
				boost::write_graphviz(os, graph, property_writer, property_writer, property_writer, vertex_id_mapper);

				lb::log_time(std::cerr) << "  Finishing output…\n";
			}

			// Sort the k-mers topologically to get the processing order.
			lb::log_time(std::cerr) << "  Sorting graph nodes topologically…\n";

			sorted_gapmers.reserve(graph.node_count());
			boost::topological_sort(
				graph,
				std::back_inserter(sorted_gapmers),
				boost::color_map(h1_graph_colour_map{graph})
			);
			std::reverse(sorted_gapmers.begin(), sorted_gapmers.end());

			// Handle the nodes in topological order.
			// semi_global_maximum_statistics will have the statistics of the semi-global maxima.
			// The algorithm works as follows:
			// – Initially colour all nodes white.
			// – Repeat until all nodes have been coloured red. (We use gapmer’s mark property for this.)
			//   – For every non-red node, follow the edges and atomically mark the target nodes grey.
			//   – For every white node in topological order, do the following:
			//     – Iterate over the edge targets and atomically update the statistics.
			//     – Mark the current (source) node red.
			//   – Colour every non-red node white.
			{
				auto const sorted_count{sorted_gapmers.size()};
				std::size_t handled_count{};
				lb::log_time(std::cerr) << "  Processing " << sorted_count << " sorted graph nodes…\n";

				// Pre-insert everything so that insertions need not be made during iteration.
				for (auto const gg : sorted_gapmers)
					semi_local_maximum_statistics.try_emplace(gg, huddinge_environment_statistics{});

				// Colour all nodes white.
#pragma omp parallel for
				for (std::size_t ii = 0; ii < sorted_count; ++ii)
				{
					auto const gg{sorted_gapmers[ii]};
					auto &count{graph.find_count(gg)};
					count.colour_white();
				}

				// Determine the maxima.
				{
					found_semi_local_maxima.clear();

#pragma omp parallel for reduction(+: handled_count)
					for (std::size_t ii = 0; ii < sorted_count; ++ii)
					{
						auto &gg{sorted_gapmers[ii]};
						auto const edge_range{graph.out_edges(gg)};
						if (edge_range.first == edge_range.second)
						{
							// The order is important b.c. mark() changes the user info bit.
							libbio_assert(not gg.is_marked());

							{
								std::lock_guard const lock{found_semi_local_maxima_mutex};
								found_semi_local_maxima.emplace_back(gg);
							}

							gg.mark();
							++handled_count;
						}
					}

					std::sort(found_semi_local_maxima.begin(), found_semi_local_maxima.end());
				}

				std::size_t round{};
				while (sorted_count - handled_count)
				{
					++round;
					if (0 == round % 1000)
						lb::log_time(std::cerr) << "   Round " << round << "; " << (sorted_count - handled_count) << " nodes remaining…\n";

					// Colour the edge targets grey.
#pragma omp parallel for schedule(static, 1)
					for (std::size_t ii = 0; ii < sorted_count; ++ii)
					{
						auto const src{sorted_gapmers[ii]};

						if (not src.is_marked())
						{
							auto edge_range{graph.out_edges(src)};
							auto edge_it{edge_range.first};
							auto const edge_end{edge_range.second};

							for (; edge_it != edge_end; ++edge_it)
							{
								auto const dst{edge_it->second};
								libbio_assert(not dst.is_marked());
								auto &dst_count{graph.find_count(dst)};
								dst_count.add_grey_colour(std::memory_order_relaxed);
							}
						}
					}

#pragma omp parallel for schedule(static, 1) reduction(+: handled_count)
					for (std::size_t ii = 0; ii < sorted_count; ++ii)
					{
						auto &src{sorted_gapmers[ii]};
						if (src.is_marked())
							continue;

						auto &src_count{graph.find_count(src)};
						if (src_count.has_grey_colour(std::memory_order_relaxed))
							continue;

						// Must be white.
						{
							auto src_stat_it{semi_local_maximum_statistics.find(src)};
							libbio_assert_neq(src_stat_it, semi_local_maximum_statistics.end());
							auto const path_length{1 + src_stat_it->second.path_length()};

							auto edge_range{graph.out_edges(src)};
							auto edge_it{edge_range.first};
							auto const edge_end{edge_range.second};

							for (; edge_it != edge_end; ++edge_it)
							{
								auto const dst{edge_it->second};
								auto dst_stat_it{semi_local_maximum_statistics.find(dst)};
								libbio_assert_neq(dst_stat_it, semi_local_maximum_statistics.end());

								// The following needs to be atomic b.c. multiple threads can access
								// the same destination at a time.
								dst_stat_it->second.update_path_length(path_length);
							}
						}

						src.mark();
						++handled_count;
					}

#pragma omp parallel for
					for (std::size_t ii = 0; ii < sorted_count; ++ii)
					{
						auto const src{sorted_gapmers[ii]};
						if (not src.is_marked())
						{
							auto &count{graph.find_count(src)};
							count.clear_grey_colour(std::memory_order_relaxed);
						}
					}
				}

				// Update the local maxima.
				{
					gapmer_huddinge_environment_statistic_map updated_statistics;
					path_length_type max_path_length{};
					for (auto const gg : found_semi_local_maxima)
					{
						auto it{semi_local_maximum_statistics.find(gg)};
						libbio_assert_neq(it, semi_local_maximum_statistics.end());

						max_path_length = std::max(max_path_length, it->second.path_length());

						auto nh{semi_local_maximum_statistics.extract(it)};
						updated_statistics.insert(std::move(nh));
					}

					using std::swap;
					swap(updated_statistics, semi_local_maximum_statistics);

					lb::log_time(std::cerr) << "   Maximum path length was " << max_path_length << ".\n";
				}
			}

			// At this point we have the local and semi-local maxima up to k - 1.
			// Update the remaining accumulators; since we only consider the values for k - 1,
			// we do not process any k-mer twice.
			auto const find_counts_and_calculate_statistics{[&](gapmer_type const gg, auto &&cb){
				auto const length{gg.length()};
				if (length == kk - 1)
				{
					auto const it{prev_counts.find(gg)};
					libbio_assert_neq(it, prev_counts.end());
					auto const counts{it->second};

					auto statistics{calculate_seed_statistics(gg, counts, kk - 1)};
					cb(statistics);
				}
			}};

			lb::log_time(std::cerr) << "  Handling local maxima…\n";
			for (auto const gg : local_maxima)
			{
				find_counts_and_calculate_statistics(gg, [&](seed_statistics const &statistics){
					accumulators.local_maximum.init_or_update(statistics);
				});
			};

			lb::log_time(std::cerr) << "  Handling semi-local maxima…\n";
			for (auto const &kv : semi_local_maximum_statistics)
			{
				auto const gg{kv.first};
				find_counts_and_calculate_statistics(gg, [&](seed_statistics &statistics){
					// Currently we only consider path length w.r.t. the semi-local Huddinge maxima.
					auto const &huddinge_env_statistics{kv.second};
					statistics.path_length = huddinge_env_statistics.path_length();

					accumulators.semi_local_maximum.init_or_update(statistics);
				});
			}

			// Report intermediate results.
			accumulators.output(std::cout, kk);

			std::cout << std::flush;
		}
	}

	return 0;
}
