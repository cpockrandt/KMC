#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include "../kmc_api/kmc_file.h"

#include <seqan3/argument_parser/all.hpp>

// #include "common.hpp"

using namespace seqan3;

int _tmain(int argc, char* argv[])
{
    argument_parser myparser{"Haplotype binning", argc, argv};

    myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Comparing databases with haplotype specific k-mers of different lengths.";
    myparser.info.version = "0.0.1";

    // int threads{omp_get_max_threads()};
    std::filesystem::path path_db_short{}, path_db_long{};

    myparser.add_positional_option(path_db_short, "Please provide the KMC db file to the shorter k-mers.");
    myparser.add_positional_option(path_db_long , "Please provide the KMC db file to the longer k-mers.");

    try
    {
        myparser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    // fundamental idea: iterate over long k-mers, and check whether it has unique short k-mers

    CKMCFile kmc_db_short;
    if (!kmc_db_short.OpenForRA(path_db_short))
    {
        std::cerr << "Could not open KMC file with short k-mers.\n";
        return 1;
    }

    CKMCFile kmc_db_long;
    if (!kmc_db_long.OpenForListing(path_db_long))
    {
        std::cerr << "Could not open KMC file with long k-mers.\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(2);

    uint32 kmer_long_length;
	{
		uint32 kmer_length;
		uint64 total_kmers;
		uint32 _mode;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint64 _max_count;
		kmc_db_short.Info(kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, total_kmers);
		std::cout << kmer_length << "-mer db size: " << total_kmers << std::endl;
		kmc_db_long.Info(kmer_long_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, total_kmers);
		std::cout << kmer_long_length << "-mer db size: " << total_kmers << std::endl;
	}

    uint64_t total_kmers = 0;
    uint64_t stats_ambiguous = 0;
    uint64_t stats_new_kmers = 0;
    uint64_t stats_contradiction = 0;

    uint32 kmer_long_c;
	// std::vector<uint64> kmer_long;
    std::vector<uint32_t> counter;
	CKmerAPI kmer_object(kmer_long_length);
	while (kmc_db_long.ReadNextKmer(kmer_object, kmer_long_c))
	{
        // TODO: can there be a short k-mer unique to mom and a short k-mer unique to dad?
        std::string const kmer_long_string{kmer_object.to_string()};
        kmc_db_short.GetCountersForRead(kmer_long_string, counter);
        // debug_stream << kmer_long_string << '\t' << kmer_long_c << '\t' << counter << '\n';

        std::array<uint32_t, 3> counter2{};
        for (uint32_t i = 0; i < counter.size(); ++i)
            ++counter2[counter[i]];

        if (counter2[0] == counter.size())
        {
            ++stats_new_kmers;
        }
        else if ((kmer_long_c == 1 && counter2[2] > 0) ||
                 (kmer_long_c == 2 && counter2[1] > 0))
        {
            ++stats_contradiction;
        }
        else if (counter2[1] > 0 && counter2[2] > 0)
        {
            ++stats_ambiguous;
        }
        // if (std::all_of(counter.begin(), counter.end(), [](uint32_t c) {return c == 0;}))
        //     ++stats_new_kmers;
        //
        // if (std::any_of(counter.begin(), counter.end(), [](uint32_t c) {return c == 1;}) &&
        //     std::any_of(counter.begin(), counter.end(), [](uint32_t c) {return c == 2;}))
        // {
        //     ++stats_ambiguous;
        // }

        ++total_kmers;

        if (total_kmers % 4096 == 0)
        {
            if (total_kmers > 4096)
                std::cout << "\033[5A";

            std::cout << "\nTotal        :\t" << total_kmers << '\n';
            std::cout << "New          :\t" << stats_new_kmers << " (" << (100.0f * stats_new_kmers / total_kmers) << " %)\n";
            std::cout << "Ambiguous    :\t" << stats_ambiguous << " (" << (100.0f * stats_ambiguous / total_kmers) << " %)\n";
            std::cout << "Contradiction:\t" << stats_contradiction << " (" << (100.0f * stats_contradiction / total_kmers) << " %)\n";
        }

		// kmer_object.to_long(kmer_long);
		// kmer_list.push_back({kmer_long[0], c});
		// if ((kmers_processed & 0b11111111111111111111) == 0)
		// {
		// 	std::cout << "\rProgress: " << (truncf((kmers_processed / kmers_total)*10000)/100) << "%   " << std::flush;
		// }
	}
	kmc_db_long.Close();

    std::cout << '\n';
    return 0;
}
