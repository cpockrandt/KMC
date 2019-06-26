#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include "../kmc_api/kmc_file.h"

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

// Get Histogram
// if (!kmer_data_base1.OpenForListing(input_file_name1))
// {
// 	print_info();
// 	return EXIT_FAILURE ;
// }
// else
// {
// 	uint32 _kmer_length1;
// 	uint32 _mode1;
//
// 	{
// 		uint32 _counter_size;
// 		uint32 _lut_prefix_length;
// 		uint32 _signature_len;
// 		uint32 _min_count;
// 		uint64 _max_count;
// 		uint64 _total_kmers;
// 		kmer_data_base1.Info(_kmer_length1, _mode1, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
// 	}
// 	std::cout << _kmer_length1 << "-mers\n";
//
// 	if(min_count_to_set)
// 	if (!(kmer_data_base1.SetMinCount(min_count_to_set)))
// 			return EXIT_FAILURE;
// 	if(max_count_to_set)
// 	if (!(kmer_data_base1.SetMaxCount(max_count_to_set)))
// 			return EXIT_FAILURE;
//
// 	uint32 c1;
// 	CKmerAPI kmer_object1(_kmer_length1);
// 	std::string kmer1;
//
// 	uint64_t hist[10000] = {0};
//
// 	while (kmer_data_base1.ReadNextKmer(kmer_object1, c1))
// 	{
// 		kmer_object1.to_string(kmer1);
//
// 		if (c1 >= 10000)
// 			c1 = 9999;
//
// 		hist[c1]++;
// 	}
//
// 	for (uint32_t i = 0; i < 10000; ++i)
// 		std::cout << i << '\t' << hist[i] << '\n';
//
// 	// fclose(out_file);
// 	kmer_data_base1.Close();
// }

int _tmain(int argc, char* argv[])
{
	argument_parser myparser{"Haplotype resolution - Finding unique k-mers", argc, argv};

	myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Determining set of unique k-mers in father and mother.";
    myparser.info.version = "0.0.1";

	uint32_t threads = omp_get_max_threads();
	std::string kmc_father_path, kmc_mother_path;

    myparser.add_positional_option(kmc_father_path, "Please provide the KMC file of the father.");
    myparser.add_positional_option(kmc_mother_path, "Please provide the KMC file of the mother.");
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

	CKMCFile kmc_db_father, kmc_db_mother;

	if (!kmc_db_father.OpenForListing(kmc_father_path))
	{
		std::cerr << "Could not open KMC file of the father.\n";
		return 1;
	}
	else if (!kmc_db_mother.OpenForListing(kmc_mother_path))
	{
		std::cerr << "Could not open KMC file of the mother.\n";
		return 1;
	}

	return 0;
}
