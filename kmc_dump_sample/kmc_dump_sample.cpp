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
	argument_parser myparser{"Haplotype resolution", argc, argv};

	myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Haplotype resolution of Trios using k-mer counting.";
    myparser.info.version = "0.0.1";

	uint32_t threads = omp_get_max_threads();
	std::string kmc_father_path, kmc_mother_path, reads_child_path;

    myparser.add_positional_option(kmc_father_path, "Please provide the KMC file of the father.");
    myparser.add_positional_option(kmc_mother_path, "Please provide the KMC file of the mother.");
    myparser.add_positional_option(reads_child_path, "Please provide the read file of the child.");
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
	// uint64_t no_reads = 0;
	uint64_t fo_mo[251-15+1][251-15+1] = {0};

	if (!kmc_db_father.OpenForRA(kmc_father_path))
	{
		std::cerr << "Could not open KMC file of the father.\n";
		return 1;
	}
	else if (!kmc_db_mother.OpenForRA(kmc_mother_path))
	{
		std::cerr << "Could not open KMC file of the mother.\n";
		return 1;
	}

	sequence_file_input fin{reads_child_path};
	auto chunks = fin
				| std::view::transform([] (auto record) { return get<field::SEQ>(record); } )
				| ranges::view::chunk(100000/1000);
	for (auto it = chunks.begin(); it != chunks.end(); ++it)
	{
		auto chunk = *it; // ranges::view::transform(*it, [](auto const & vec) { return ranges::view::transform(vec, [](dna5 c) { return c.to_char(); } ); });
		std::vector<std::string> chunks_container;
		chunks_container.reserve(100000/1000);
		for (auto chunk_it = chunk.begin(); chunk_it != chunk.end(); chunk_it++)
		{
			auto vec = *chunk_it;
			if (vec.size() == 250)
			{
				auto read = vec | view::to_char;
				chunks_container.emplace_back(std::string(read.begin(), read.end()));
			}
		}

		// #pragma omp parallel for num_threads(threads)
		for (auto x = chunks_container.begin(); x < chunks_container.end(); ++x)
		{
			std::vector<uint32> counters_father, counters_mother;
			kmc_db_father.GetCountersForRead_kmc2(*x, counters_father);
			debug_stream << counters_father << std::endl;
			kmc_db_mother.GetCountersForRead_kmc2(*x, counters_mother);
			debug_stream << counters_mother << std::endl;

			uint32_t mo = 0, fo = 0;
			for (uint32_t i = 0; i < counters_father.size(); ++i)
			{
				if (counters_father[i] < 5 && counters_mother[i] >= 5)
					++mo;
				else if (counters_father[i] >= 5 && counters_mother[i] < 5)
					++fo;
			}

			std::cout << mo << ' ' << fo << std::endl << std::endl;

			#pragma omp atomic
			fo_mo[fo][mo]++;
			// #pragma omp atomic
			// ++no_reads;

			// if (no_reads > 5)
			// 	return 0;

			// if ((no_reads & 0b11111111111111111) == 0b11111111111111111)
			// 	std::cout << '.' << std::flush;
		}

	}

	std::cout << '\n' << no_reads << '\n';

	for (uint32_t i = 0; i < 251; ++i)
	{
		for (uint32_t j = 0; j < 251; ++j)
		{
			std::cout << fo_mo[i][j] << '\t';
		}
		std::cout << '\n';
	}

	return 0;

	// if (!kmer_data_base1.OpenForRA(input_file_name1))
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
	// 	kmer1 = "AACG";
	// 	std::cout << kmer1 << '\n';
	// 	kmer_object1.from_string(kmer1.c_str());
	// 	if (kmer_data_base1.CheckKmer(kmer_object1, c1))
	// 		std::cout << "Found: " << c1 << '\n';
	// 	else
	// 		std::cout << "Not found!\n";
	//
	// 	kmer_data_base1.Close();
	// }
}
