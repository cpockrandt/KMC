#include "omp.h"
#include <algorithm>
#include <iostream>
#include <set>
#include <inttypes.h>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/range/view/all.hpp>

using namespace seqan3;

bool checkRead(auto const & read, auto const & found_kmers, std::vector<uint64_t> & new_kmers)
{
	bool hit = false;

	auto my_rank = [](dna5 const base)
	{
		return static_cast<dna4>(base);
		// switch (base)
		// {
		// 	case
		// }
	};

	auto const & ranks  = get<field::SEQ>(read)
					    | ranges::view::transform(my_rank)
						| seqan3::view::to_rank;
	auto const & ranksC = get<field::SEQ>(read)
						| ranges::view::transform(my_rank)
						| seqan3::view::complement
						| seqan3::view::to_rank; // this is just the complement, not the reverse complement!

 	// debug_stream << "Seq    : " << get<field::SEQ>(read) << std::endl;
	// debug_stream << "Ranks  : " << ranks << std::endl;
	// debug_stream << "RanksC : " << ranksC << std::endl;

	if (ranks.size() >= 31)
	{
		uint64_t kmer = 0, kmerRC = 0;
		uint32_t i = 0;
		for (; i < 31; ++i)
		{
			kmer <<= 2;
			kmer |= ranks[i];
			kmerRC |= (static_cast<uint64_t>(ranksC[i]) << (2 * i)); // fill it the other way around to reverse the complementary seq.
		}
		uint64_t const canonical_kmer = std::min(kmer, kmerRC);
		// debug_stream << "kmer   : " << std::bitset<64>(kmer) << std::endl;
		// debug_stream << "kmerRC : " << std::bitset<64>(kmerRC) << std::endl;
		if (found_kmers.find(canonical_kmer) == found_kmers.end())
		{
			new_kmers.push_back(canonical_kmer);
		}
		else
		{
			// debug_stream << "Hit1!" << std::endl;
			hit = true;
		}

		for (; i < ranks.size(); ++i)
		{
			kmer <<= 2;
			kmer &= (static_cast<uint64_t>(1) << 62) - 1; // 2^62 - 1, k-mer only uses rightmost 2*31 bits
			kmer |= ranks[i];
			kmerRC >>= 2;
			kmerRC |= (static_cast<uint64_t>(ranksC[i]) << (2 * 30));
			uint64_t const canonical_kmer = std::min(kmer, kmerRC);
			if (found_kmers.find(canonical_kmer) == found_kmers.end())
			{
				new_kmers.push_back(canonical_kmer);
			}
			else
			{
				// debug_stream << "Hit2!" << std::endl;
				hit = true;
			}
			// debug_stream << "kmer   : " << std::bitset<64>(kmer) << std::endl;
			// debug_stream << "kmerRC : " << std::bitset<64>(kmerRC) << std::endl;
		}
	}

	return hit;
}

int main(int argc, char* argv[])
{
	argument_parser myparser{"Comparing read data sets.", argc, argv};

	myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Comparing read data sets.";
    myparser.info.version = "0.0.1";

	std::string reads_path1, reads_path2;
	uint32_t threads = omp_get_max_threads();

    myparser.add_positional_option(reads_path1, "Please provide a read data set (PE 1).");
    myparser.add_positional_option(reads_path2, "Please provide a read data set (PE 2).");
    myparser.add_option(threads, 't', "threads", "Number of threads", option_spec::DEFAULT, arithmetic_range_validator{1, 64});

    try
    {
         myparser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

	sequence_file_input reads1{reads_path1}, reads2{reads_path2};
	std::vector<typename seqan3::sequence_file_input<>::record_type> reads1_container, reads2_container;

	#pragma omp parallel num_threads(2)
	{
		#pragma omp sections
		{
			#pragma omp section
			{
				reads1_container.reserve(3'000'000);
				for (auto it = reads1.begin(); it != reads1.end(); ++it)
					reads1_container.push_back(*it);
				std::cout << reads1_container.size() << " reads read (PE 1)." << std::endl;
			}

			#pragma omp section
			{
				reads2_container.reserve(3'000'000);
				for (auto it = reads2.begin(); it != reads2.end(); ++it)
					reads2_container.push_back(*it);
				std::cout << reads2_container.size() << " reads read (PE 2)." << std::endl;
			}
		}
	}

	if (reads1_container.size() != reads2_container.size())
	{
		std::cerr << "Different number of PE reads ...\n";
		exit(1);
	}

	std::set<uint64_t> found_kmers;

	uint64_t readno = 0;
	uint64_t someNreadno = 0;
	uint64_t allNreadno = 0;
	uint64_t uniqueKreadno = 0;

    sequence_file_output fout_unique1{"unique.1.fq"}, fout_unique2{"unique.2.fq"};

	// #pragma omp parallel for num_threads(threads)
	for (uint64_t i = 0; i < reads1_container.size(); ++i)
	{
		auto const & read1 = get<field::SEQ>(reads1_container[i]);
		auto const & read2 = get<field::SEQ>(reads2_container[i]);

		// #pragma omp atomic
		++readno;

		auto isN = [](auto const & base) { return base == 'N'_dna5; };

		// remove reads that only consist of Ns
		if (std::all_of(read1.begin(), read1.end(), isN) && std::all_of(read2.begin(), read2.end(), isN))
		{
			// #pragma omp atomic
			++allNreadno;
		}
		// keep reads that have at least one N
		else if (std::any_of(read1.begin(), read1.end(), isN) || std::any_of(read2.begin(), read2.end(), isN))
		{
			// #pragma omp atomic
			++someNreadno;

			// #pragma omp critical
			{
				fout_unique1.push_back(reads1_container[i]);
				fout_unique2.push_back(reads2_container[i]);
			}
		}
		else
		{
			std::vector<uint64_t> new_kmers; // TODO: move outside of loop
			bool hit = checkRead(reads1_container[i], found_kmers, new_kmers);
			if (checkRead(reads2_container[i], found_kmers, new_kmers))
			{
				hit = true;
			}

			if (!hit)
			{
				// #pragma omp critical // TODO: keep an output buffer for each thread
				{
					++uniqueKreadno;
					fout_unique1.push_back(reads1_container[i]);
					fout_unique2.push_back(reads2_container[i]);
				}
			}

			// #pragma omp critical
			{
				for (auto const can_kmer : new_kmers)
					found_kmers.insert(can_kmer);
			}
		}


		if ((readno & 0b11111111111111111) == 0)
		{
			std::cout << "\rProgress: " << (truncf((readno / (reads1_container.size()*1.0f))*10000)/100) << "%   " << std::flush;
		}
	}

	std::cout << "pure N reads: " << allNreadno << " (removed)\n";
	std::cout << "some N reads: " << someNreadno << " (kept)\n";
	std::cout << "unique k-reads: " << uniqueKreadno << " (kept)\n";

	return 0;
}
