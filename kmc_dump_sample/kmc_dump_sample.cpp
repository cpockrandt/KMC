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

	uint32_t min = 3;
	uint32_t threads = omp_get_max_threads();
	std::filesystem::path kmc_father_path{}, kmc_mother_path{}, reads_child_path1{}, reads_child_path2{}, out_dir{};

    myparser.add_positional_option(kmc_father_path, "Please provide the KMC file of the father.");
    myparser.add_positional_option(kmc_mother_path, "Please provide the KMC file of the mother.");
    myparser.add_positional_option(reads_child_path1, "Please provide the read file of the child (1. PE).");
    myparser.add_positional_option(reads_child_path2, "Please provide the read file of the child (2. PE).");
    myparser.add_option(threads, 't', "threads", "Number of threads", option_spec::DEFAULT, arithmetic_range_validator{1, 64});
    myparser.add_option(min, 'm', "min", "Min k-mers (occurring less than X times)", option_spec::DEFAULT);
	myparser.add_option(out_dir, 'o', "out", "Directory to output binned reads.", option_spec::DEFAULT, output_directory_validator{});

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

	uint32 kmer_length;
	{
		uint32 _mode;
		uint32 _kmer_length2;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint64 _max_count;
		uint64 _total_kmers;
		kmc_db_father.Info(kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
		kmc_db_mother.Info(_kmer_length2, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

		if (kmer_length != _kmer_length2)
		{
			std::cerr << "The databases have different k-mer lengths.\n";
			return 1;
		}
		std::cout << "k-mer length: " << kmer_length << std::endl;
	}

	uint64_t no_reads = 0;
	uint64_t fo_mo[2 * (250 - kmer_length + 1) + 1][2 * (250 - kmer_length + 1) + 1] = {0};

	sequence_file_input fin1{reads_child_path1}, fin2{reads_child_path2};
	auto chunks1 = fin1
				 | std::view::transform([] (auto record) { return get<field::SEQ>(record); } )
				 | ranges::view::chunk(500000);
	auto chunks2 = fin2
				 | std::view::transform([] (auto record) { return get<field::SEQ>(record); } )
				 | ranges::view::chunk(500000);

	std::vector</*std::string*/ typename sequence_file_input<>::record_type> chunks_container1, chunks_container2;
	chunks_container1.reserve(500000);
	chunks_container2.reserve(500000);

	// TODO: automatic file ending
    sequence_file_output fout_mother_only1{out_dir / "mother_only.1.fa"}, fout_mother_only2{out_dir / "mother_only.2.fa"},
						 fout_father_only1{out_dir / "father_only.1.fa"}, fout_father_only2{out_dir / "father_only.2.fa"},
						 fout_mixed1{out_dir / "mixed.1.fa"},             fout_mixed2{out_dir / "mixed.2.fa"};

	auto it1{chunks1.begin()};
	auto it2{chunks2.begin()};
	for (; it1 != chunks1.end(); ++it1, ++it2)
	{
		// auto read = vec | view::to_char;
		// chunks_container.emplace_back(std::string(read.begin(), read.end()));
		#pragma omp parallel num_threads(2)
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					chunks_container1.clear();
					auto chunk = *it1; // ranges::view::transform(*it, [](auto const & vec) { return ranges::view::transform(vec, [](dna5 c) { return c.to_char(); } ); });
					for (auto chunk_it = chunk.begin(); chunk_it != chunk.end(); chunk_it++)
					{
						auto vec = *chunk_it;
						chunks_container1.push_back(vec);
					}
				}

				#pragma omp section
				{
					chunks_container2.clear();
					auto chunk = *it2; // ranges::view::transform(*it, [](auto const & vec) { return ranges::view::transform(vec, [](dna5 c) { return c.to_char(); } ); });
					for (auto chunk_it = chunk.begin(); chunk_it != chunk.end(); chunk_it++)
					{
						auto vec = *chunk_it;
						chunks_container2.emplace_back(*chunk_it);
					}
				}
			}
		}

		#pragma omp parallel for num_threads(threads)
		for (uint64_t id = 0; id < chunks_container1.size(); ++id)
		{
			auto const & read1 = chunks_container1[id];
			auto const & read2 = chunks_container2[id];

			auto read1_view = get<field::SEQ>(read1) | view::to_char;
			auto read2_view = get<field::SEQ>(read2) | view::to_char;
			std::string read1_str(read1_view.begin(), read1_view.end());
			std::string read2_str(read2_view.begin(), read2_view.end());

			std::vector<uint32> counters_father1, counters_mother1, counters_father2, counters_mother2; // TODO: put outside of loop
			kmc_db_father.GetCountersForRead(read1_str, counters_father1);
			kmc_db_mother.GetCountersForRead(read1_str, counters_mother1);
			kmc_db_father.GetCountersForRead(read2_str, counters_father2);
			kmc_db_mother.GetCountersForRead(read2_str, counters_mother2);
			// debug_stream << counters_father << std::endl;
			// debug_stream << counters_mother << std::endl;

			uint32_t mo = 0, fo = 0;
			for (uint32_t i = 0; i < counters_father1.size(); ++i)
			{
				if (counters_father1[i] < min && counters_mother1[i] >= min)
					++mo;
				else if (counters_father1[i] >= min && counters_mother1[i] < min)
					++fo;
			}
			for (uint32_t i = 0; i < counters_father2.size(); ++i)
			{
				if (counters_father2[i] < min && counters_mother2[i] >= min)
					++mo;
				else if (counters_father2[i] >= min && counters_mother2[i] < min)
					++fo;
			}
			// std::cout << mo << ' ' << fo << std::endl << std::endl;

			if (mo > 0 && fo == 0)
			{
				#pragma omp critical(fout_mother_only)
				{
					fout_mother_only1.push_back(read1);
					fout_mother_only2.push_back(read2);
				}
			}
			else if (fo > 0 && mo == 0)
			{
				#pragma omp critical(fout_father_only)
				{
					fout_father_only1.push_back(read1);
					fout_father_only2.push_back(read2);
				}
			}
			else
			{
				#pragma omp critical(fout_mixed)
				{
					fout_mixed1.push_back(read1);
					fout_mixed2.push_back(read2);
				}
			}

			#pragma omp atomic
			fo_mo[fo][mo]++;
			#pragma omp atomic
			++no_reads;

			// if ((no_reads & 0b11111111111111111) == 0b11111111111111111)
			// 	std::cout << '.' << std::flush;
		}
		std::cout << '.' << std::flush;

	}

	// std::cout << '\n' << no_reads << '\n';

	// determine last row that contains non-zero values
	uint32_t last_row = 0;
	for (uint32_t j = 0; j < 2 * (250 - kmer_length + 1) + 1; ++j)
	{
		for (uint32_t i = 0; i < 2 * (250 - kmer_length + 1) + 1; ++i)
		{
			if (fo_mo[i][j] > 0 && i > last_row)
			{
				last_row = i;
			}
		}
	}
	std::cout << "Last row: " << last_row << '\n';

	for (uint32_t i = 0; i < last_row + 1; ++i)
	{
		uint32_t last_col = 0;
		for (uint32_t j = 0; j < 2 * (250 - kmer_length + 1) + 1; ++j)
		{
			if (fo_mo[i][j] > 0 && j > last_col)
			{
				last_col = j;
			}
		}

		for (uint32_t j = 0; j < last_col + 1; ++j)
		{
			if (fo_mo[i][j] > 0)
				std::cout << fo_mo[i][j] << '\t';
			else
				std::cout << '\t';
		}
		std::cout << '\n';
	}

	uint64_t fa_only = 0;
	uint64_t mo_only = 0;
	uint64_t mixed = 0;
	for (uint32_t i = 1; i < 2 * (250 - kmer_length + 1) + 1; ++i)
	{
		fa_only += fo_mo[i][0];
		mo_only += fo_mo[0][i];

		for (uint32_t j = 1; j < 2 * (250 - kmer_length + 1) + 1; ++j)
		{
			mixed += fo_mo[i][j];
		}
	}

	std::cout << "\n\nTotal:\t" << no_reads << '\n';
	std::cout << "None:\t" << fo_mo[0][0] << '\n';
	std::cout << "Father only:\t" << fa_only << '\n';
	std::cout << "Mother only:\t" << mo_only << '\n';
	std::cout << "Mixed:\t" << mixed << '\n';

	return 0;
}
