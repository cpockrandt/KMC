/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  This file demonstrates the example usage of kmc_api software.
  It reads kmer_counter's output and prints kmers to an output file.

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.1.1
  Date   : 2019-05-19
*/

#include "stdafx.h"
#include <iostream>
#include "../kmc_api/kmc_file.h"

void print_info(void);

//----------------------------------------------------------------------------------
// Check if --help or --version was used
bool help_or_version(int argc, char** argv)
{
	const std::string version = "--version";
	const std::string help = "--help";
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == version || argv[i] == help)
			return true;
	}
	return false;
}

int _tmain(int argc, char* argv[])
{
	if (argc == 1 || help_or_version(argc, argv))
	{
		print_info();
		return 0;
	}

	CKMCFile kmer_data_base1, kmer_data_base2;
	int32 i;
	uint32 min_count_to_set = 0;
	uint32 max_count_to_set = 0;
	std::string input_file_name1, input_file_name2;
	std::string output_file_name;

	// FILE * out_file;
	//------------------------------------------------------------
	// Parse input parameters
	//------------------------------------------------------------
	if(argc < 3)
	{
		print_info();
		return EXIT_FAILURE;
	}

	for(i = 1; i < argc; ++i)
	{
		if(argv[i][0] == '-')
		{
			if(strncmp(argv[i], "-ci", 3) == 0)
				min_count_to_set = atoi(&argv[i][3]);
			else if(strncmp(argv[i], "-cx", 3) == 0)
					max_count_to_set = atoi(&argv[i][3]);
		}
		else
			break;
	}

	if(argc - i < 2)
	{
		print_info();
		return EXIT_FAILURE;
	}

	input_file_name1 = std::string(argv[i++]);
	output_file_name = std::string(argv[i]);

	// if((out_file = fopen (output_file_name.c_str(),"wb")) == NULL)
	// {
	// 	print_info();
	// 	return EXIT_FAILURE;
	// }
	//
	// setvbuf(out_file, NULL ,_IOFBF, 1 << 24);

	//------------------------------------------------------------------------------
	// Open kmer database for listing and print kmers within min_count and max_count
	//------------------------------------------------------------------------------
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

	// TODO: read fastq.gz file with SeqAn3
	// std::vector<uint32> counters1, counters2;
	// kmer_data_base1.GetCountersForRead_kmc2(std::string read, counters1);
	// kmer_data_base2.GetCountersForRead_kmc2(std::string read, counters2);
	uint32_t noMatches1 = std::count_if(counters1.begin(), counters1.end(), [](uint32 i){ return i < 5 && i < 45; });
	uint32_t noMatches2 = std::count_if(counters2.begin(), counters2.end(), [](uint32 i){ return i < 5 && i < 45; });

	if (!kmer_data_base1.OpenForRA(input_file_name1))
	{
		print_info();
		return EXIT_FAILURE ;
	}
	else
	{
		uint32 _kmer_length1;
		uint32 _mode1;

		{
			uint32 _counter_size;
			uint32 _lut_prefix_length;
			uint32 _signature_len;
			uint32 _min_count;
			uint64 _max_count;
			uint64 _total_kmers;
			kmer_data_base1.Info(_kmer_length1, _mode1, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
		}
		std::cout << _kmer_length1 << "-mers\n";

		if(min_count_to_set)
		if (!(kmer_data_base1.SetMinCount(min_count_to_set)))
				return EXIT_FAILURE;
		if(max_count_to_set)
		if (!(kmer_data_base1.SetMaxCount(max_count_to_set)))
				return EXIT_FAILURE;

		uint32 c1;
		CKmerAPI kmer_object1(_kmer_length1);
		std::string kmer1;

		kmer1 = "AACG";
		std::cout << kmer1 << '\n';
		kmer_object1.from_string(kmer1.c_str());
		if (kmer_data_base1.CheckKmer(kmer_object1, c1))
			std::cout << "Found: " << c1 << '\n';
		else
			std::cout << "Not found!\n";

		kmer_data_base1.Close();
	}

	// bool CKMCFile::GetCountersForRead_kmc2(const std::string& read, std::vector<uint32>& counters)

	return EXIT_SUCCESS;
}
// -------------------------------------------------------------------------
// Print execution options
// -------------------------------------------------------------------------
void print_info(void)
{
	std::cout << "KMC dump ver. " << KMC_VER << " (" << KMC_DATE << ")\n"
			  << "\nUsage:\nkmc_dump [options] <kmc_database> <output_file>\n"
			  << "Parameters:\n"
			  << "<kmc_database> - kmer_counter's output\n"
			  << "Options:\n"
			  << "-ci<value> - exclude k-mers occurring less than <value> times\n"
			  << "-cx<value> - exclude k-mers occurring more of than <value> times\n";
};

// ***** EOF
