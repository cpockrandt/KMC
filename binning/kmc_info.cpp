#include "omp.h"
#include "../kmc_dump_sample/stdafx.h"
#include <iostream>
#include <iomanip>
#include "../kmc_api/kmc_file.h"

#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

int _tmain(int argc, char* argv[])
{
    argument_parser parser{"KMC-info", argc, argv};

    parser.info.author = "Christopher Pockrandt";
    parser.info.short_description = "Get info on KMC-db.";
    parser.info.version = "0.0.1";

    std::filesystem::path path;
    parser.add_positional_option(path, "Prefix of KMC database.");

    bool merged_db;
    parser.add_flag(merged_db, 'm', "merged", "Merged database (father and mother)");

    try
    {
        parser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    uint32_t kmer_length;
    uint64 unique_kmers;

    {
	CKMCFile kmc;
        if (!kmc.OpenForRA(path))
	{
	    std::cerr << "Could not open KMC file for RA.\n";
	    return 1;
	}
        uint32 _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
        uint64 _max_count;
        kmc.Info(kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, unique_kmers);
        std::cout << "kmer-length:  " << kmer_length << '\n'
		  << "counter-size: " << _counter_size << '\n'
		  << "min count:    " << _min_count << '\n'
		  << "max count:    " << _max_count << '\n'
		  << "unique kmers: " << unique_kmers << '\n';
    }

    //std::cout << "Start counting total kmers ...\n";

    CKMCFile kmc;
    if (!kmc.OpenForListing(path))
    {
        std::cerr << "Could not open KMC for Listing.\n";
	return 1;
    }

    uint32_t c;
    uint64_t total_kmers = 0, mom = 0, dad = 0;
    CKmerAPI kmer_object(kmer_length);
    while (kmc.ReadNextKmer(kmer_object, c))
    {
	total_kmers += c;
	if (merged_db && c == 1)
	    ++dad;
        else if (merged_db && c == 2)
            ++mom;
    }

    if (merged_db)
    {
        std::cout << "Father:       " << dad << '\n';
        std::cout << "Mother:       " << mom << '\n';
    }
    else
    {
        std::cout << "Total kmers:  " << total_kmers << '\n';
    }

    return 0;
}
