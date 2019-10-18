#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include "../kmc_api/kmc_file.h"

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

using kmer_repr_t = uint64_t;

inline void read_kmc_db(auto & kmc_db, auto & kmer_list, uint32 const kmer_length) noexcept
{
	uint32 c;
	std::vector<uint64> kmer_long;
	CKmerAPI kmer_object(kmer_length);
	while (kmc_db.ReadNextKmer(kmer_object, c))
	{
		kmer_object.to_long(kmer_long);
		kmer_list.push_back(kmer_long[0]);
	}
	kmc_db.Close();
}

int _tmain(int argc, char* argv[])
{
	argument_parser myparser{"Comparing read data sets.", argc, argv};

	myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Comparing read data sets.";
    myparser.info.version = "0.0.1";

	// uint32_t min = 0;
	std::string kmc_path, reads_path1, reads_path2;
	uint32_t threads = omp_get_max_threads();

    myparser.add_positional_option(kmc_path, "Please provide the KMC file of a read dataset.");
    myparser.add_positional_option(reads_path1, "Please provide a read data set (PE 1).");
    myparser.add_positional_option(reads_path2, "Please provide a read data set (PE 2).");
    // myparser.add_option(min, 'm', "min", "Min k-mers (occurring less than X times)", option_spec::DEFAULT);
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

	CKMCFile kmc_db;

	if (!kmc_db.OpenForRA(kmc_path))
	{
		std::cerr << "Could not open KMC file of the healthy read dataset.\n";
		return 1;
	}

	{
		uint32 kmer_length;
		uint64 total_kmers;
		uint32 _mode;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint64 _max_count;
		kmc_db.Info(kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, total_kmers);

		std::cout << kmer_length << "-mer db size: " << total_kmers << std::endl;
	}

	uint64_t readCount1[121 + 1] = {0};
	uint64_t readCount2[121 + 1] = {0};

	sequence_file_input reads1{reads_path1}, reads2{reads_path2};
	// auto reads1 = fin1 | std::views::transform([] (auto record) { return get<field::SEQ>(record); } );
	// auto reads2 = fin2 | std::views::transform([] (auto record) { return get<field::SEQ>(record); } );
	std::vector<typename seqan3::sequence_file_input<>::record_type> reads1_container, reads2_container;

	#pragma omp parallel num_threads(2)
	{
		#pragma omp sections
		{
			#pragma omp section
			{
				reads1_container.reserve(80'000'000);
				for (auto it = reads1.begin(); it != reads1.end(); ++it)
					reads1_container.push_back(*it);
				std::cout << reads1_container.size() << " reads read (PE 1)." << std::endl;
			}

			#pragma omp section
			{
				reads2_container.reserve(80'000'000);
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

	std::vector<std::vector<uint32>> counters1, counters2;
	counters1.resize(threads);
	counters2.resize(threads);

	uint64_t readno = 0;
	uint64_t readno_onlyNs = 0;
	uint64_t readno_set_difference = 0;

    sequence_file_output fout1_0{"out.nohit.1.fq"}, fout2_0{"out.nohit.2.fq"};
    //sequence_file_output fout1_10{"out.fewhits.1.fq"}, fout2_10{"out.fewhits.2.fq"};

	#pragma omp parallel for num_threads(threads)
	for (uint64_t i = 0; i < reads1_container.size(); ++i)
	{
		auto read1 = get<field::SEQ>(reads1_container[i]) | views::to_char;
		auto read2 = get<field::SEQ>(reads2_container[i]) | views::to_char;

		#pragma omp atomic
		++readno;

		if (std::all_of(read1.begin(), read1.end(), [](auto const & base) { return base == 'N'; }) &&
			std::all_of(read2.begin(), read2.end(), [](auto const & base) { return base == 'N'; }))
		{
			#pragma omp atomic
			++readno_onlyNs;
		}
		else
		{
			auto & counter1 = counters1[omp_get_thread_num()];
			auto & counter2 = counters2[omp_get_thread_num()];

			kmc_db.GetCountersForRead_kmc2_both_strands(std::string(read1.begin(), read1.end()), counter1);
			kmc_db.GetCountersForRead_kmc2_both_strands(std::string(read2.begin(), read2.end()), counter2);
			// debug_stream << counters << std::endl;

			uint16_t c1 = std::min<uint16_t>(std::count_if(counter1.begin(), counter1.end(), [](auto i) { return i > 0; }), 121);
			uint16_t c2 = std::min<uint16_t>(std::count_if(counter2.begin(), counter2.end(), [](auto i) { return i > 0; }), 121);

			#pragma omp atomic
			readCount1[c1]++;
			#pragma omp atomic
			readCount2[c2]++;

			/*if (c1 < 10 || c2 < 10)
			{
				#pragma omp critical
				{
					fout1_10.push_back(reads1_container[i]);
					fout2_10.push_back(reads2_container[i]);
				}
			}*/

			if (c1 == 0 && c2 == 0)
			{
				#pragma omp critical
				{
					fout1_0.push_back(reads1_container[i]);
					fout2_0.push_back(reads2_container[i]);
					++readno_set_difference;
				}
			}
		}

		if ((readno & 0b11111111111111111) == 0)
		{
			std::cout << "\rProgress: " << (truncf((readno / (reads1_container.size()*1.0f))*10000)/100) << "%   " << std::flush;
		}
	}

	std::cout << "Reads processed : " << readno << '\n'
	 		  << "Reads eliminated: " << readno_onlyNs << " (weird NNNNN....N reads)\n"
	 		  << "Reads outputted : " << readno_set_difference << " (reads in set difference, in out.nohit.[1/2].fq)\n";

	// std::cout << "Distribution of bins (0 k-mer matches to 121 k-mer matches per read. For 1.fq and 2.fq separately." << std::endl;
	//
	// for (uint16_t i = 0; i < 121 + 1; ++i)
	// 	std::cout << readCount1[i] << ';';
	// std::cout << std::endl << std::endl;
	//
	// for (uint16_t i = 0; i < 121 + 1; ++i)
	// 	std::cout << readCount2[i] << ';';
	// std::cout << std::endl;

	return 0;
}
