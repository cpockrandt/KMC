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
	input_file_name2 = std::string(argv[i++]);
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

	if (!kmer_data_base1.OpenForListing(input_file_name1) || !kmer_data_base2.OpenForListing(input_file_name2))
	{
		print_info();
		return EXIT_FAILURE ;
	}
	else
	{
		uint32 _kmer_length1, _kmer_length2;
		uint32 _mode1, _mode2;

		{
			uint32 _counter_size;
			uint32 _lut_prefix_length;
			uint32 _signature_len;
			uint32 _min_count;
			uint64 _max_count;
			uint64 _total_kmers;
			kmer_data_base1.Info(_kmer_length1, _mode1, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
			kmer_data_base2.Info(_kmer_length2, _mode2, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

			if (_mode1 != _mode2 || _kmer_length1 != _kmer_length2)
				exit(42);
		}



		if(min_count_to_set)
		if (!(kmer_data_base1.SetMinCount(min_count_to_set)) || !(kmer_data_base2.SetMinCount(min_count_to_set)))
				return EXIT_FAILURE;
		if(max_count_to_set)
		if (!(kmer_data_base1.SetMaxCount(max_count_to_set)) || !(kmer_data_base2.SetMaxCount(max_count_to_set)))
				return EXIT_FAILURE;

		// if (_mode1) //quake compatible mode
		// {
		// 	float counter;
		// 	while (kmer_data_base.ReadNextKmer(kmer_object, counter))
		// 	{
		// 		kmer_object.to_string(str);
		// 		// fprintf(out_file, "%s\t%f\n", str.c_str(), counter);
		// 		std::cout << str.c_str() << '\t' << counter << std::endl;
		// 	}
		// }
		// else
		{
			uint32 c1, c2;
			CKmerAPI kmer_object1(_kmer_length1), kmer_object2(_kmer_length2);
			std::string kmer1, kmer2;

			bool db1_empty = !kmer_data_base1.ReadNextKmer(kmer_object1, c1);
			kmer_object1.to_string(kmer1);
			bool db2_empty = !kmer_data_base2.ReadNextKmer(kmer_object2, c2);
			kmer_object2.to_string(kmer2);

			uint64_t no_equal = 0;

			uint64_t no_not_equal = 0;
			uint64_t sum_not_equal = 0;

			uint64_t no_only_in_db1 = 0;
			uint64_t sum_only_in_db1 = 0;

			uint64_t no_only_in_db2 = 0;
			uint64_t sum_only_in_db2 = 0;

			while (!db1_empty || !db2_empty)
			{
				if (kmer1 == kmer2) // TODO: compare hashes
				{
					if (c1 == c2)
						++no_equal;
					else
					{
						++no_not_equal;
						sum_not_equal += std::max(c1, c2) - std::min(c1, c2);
					}

					// std::cout << kmer1 << '\n';
					db1_empty = !kmer_data_base1.ReadNextKmer(kmer_object1, c1);
					kmer_object1.to_string(kmer1);
					db2_empty = !kmer_data_base2.ReadNextKmer(kmer_object2, c2);
					kmer_object2.to_string(kmer2);
				}
				else if (kmer1 > kmer2)
				{
					// k-mer only occurs in DB2 ('c2' times)
					++no_only_in_db2;
					sum_only_in_db2 += c2;

					// std::cout << kmer2 << '\n';
					db2_empty = !kmer_data_base2.ReadNextKmer(kmer_object2, c2);
					kmer_object2.to_string(kmer2);
				}
				else // if (kmer1 < kmer2)
				{
					// k-mer only occurs in DB1 ('c1' times)
					++no_only_in_db1;
					sum_only_in_db1 += c1;

					// std::cout << kmer1 << '\n';
					db1_empty = !kmer_data_base1.ReadNextKmer(kmer_object1, c1);
					kmer_object1.to_string(kmer1);
				}

				if (db1_empty)
					kmer1 = std::string(_kmer_length1, 'Z');
				if (db2_empty)
					kmer2 = std::string(_kmer_length2, 'Z');
			}

			std::cout <<
				"k-mer      : " << _kmer_length1 << '\n' <<
				"equal      : " << no_equal << '\n' <<
				"not equal  : " << no_not_equal << " (" << (1.0d * sum_not_equal / no_not_equal) << ")" << '\n' <<
				"only in db1: " << no_only_in_db1 << " (" << (1.0d * sum_only_in_db1 / no_only_in_db1) << ")" << '\n' <<
				"only in db2: " << no_only_in_db2 << " (" << (1.0d * sum_only_in_db2 / no_only_in_db2) << ")" << '\n';

			// while (kmer_data_base1.ReadNextKmer(kmer_object1, counter))
			// {
			// 	std::cout << str.c_str() << '\t' << counter << std::endl;
			// 	// fprintf(out_file, "%s\t%u\n", str.c_str(), counter);
			// }
			// std::cout << "------------------------------------------------\n";
			// while (kmer_data_base2.ReadNextKmer(kmer_object2, counter))
			// {
			// 	kmer_object2.to_string(str);
			// 	std::cout << str.c_str() << '\t' << counter << std::endl;
			// 	// fprintf(out_file, "%s\t%u\n", str.c_str(), counter);
			// }
		}


		// fclose(out_file);
		kmer_data_base1.Close();
		kmer_data_base2.Close();
	}

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
