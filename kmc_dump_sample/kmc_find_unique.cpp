#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include "../kmc_api/kmc_file.h"

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

using kmer_repr_t = uint64_t;

inline void read_kmc_db(auto & kmc_db, auto & kmer_list, uint32 const kmer_length, auto & kmer_list2, double const kmers_total) noexcept
{
	uint32 c;
	std::vector<uint64> kmer_long;
	CKmerAPI kmer_object(kmer_length);
	while (kmc_db.ReadNextKmer(kmer_object, c))
	{
		kmer_object.to_long(kmer_long);
		kmer_list.push_back({kmer_long[0], c});

		uint64_t const kmers_processed = kmer_list.size() + kmer_list2.size();
		if ((kmers_processed & 0b11111111111111111111) == 0)
		{
			// #pragma omp single
			std::cout << "\rProgress: " << (truncf((kmers_processed / kmers_total)*10000)/100) << "%   " << std::flush;
		}
	}
	kmc_db.Close();
}

int _tmain(int argc, char* argv[])
{
	argument_parser myparser{"Haplotype resolution - Finding unique k-mers", argc, argv};

	myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Determining set of unique k-mers in father and mother.";
    myparser.info.version = "0.0.1";

	uint32_t min = 0;
	std::string kmc_father_path{}, kmc_mother_path{}, out_dir{};

    myparser.add_positional_option(kmc_father_path, "Please provide the KMC file of the father.");
    myparser.add_positional_option(kmc_mother_path, "Please provide the KMC file of the mother.");
    myparser.add_option(min, 'm', "min", "Min k-mers (occurring less than X times)", option_spec::DEFAULT);

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
	else if (!kmc_db_father.SetMinCount(min))
	{
		std::cerr << "Cannot set min-count on KMC file of the father.\n";
		return 1;
	}
	else if (!kmc_db_mother.SetMinCount(min))
	{
		std::cerr << "Cannot set min-count on KMC file of the mother.\n";
		return 1;
	}

	uint32 kmer_length;
	uint64 total_kmers1, total_kmers2;
	{
		uint32 _mode;
		uint32 _kmer_length2;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint64 _max_count;
		kmc_db_father.Info(kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, total_kmers1);
		kmc_db_mother.Info(_kmer_length2, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, total_kmers2);

		if (kmer_length != _kmer_length2)
		{
			std::cerr << "The databases have different k-mer lengths.\n";
			return 1;
		}
		std::cout << "k-mer db sizes (father and mother, without -m): " << total_kmers1 << " and " << total_kmers2 << std::endl;
	}

	std::vector<std::pair<kmer_repr_t, uint32_t>> kmer_list_father, kmer_list_mother;

	struct {
		bool operator()(std::pair<kmer_repr_t, uint32_t> a, std::pair<kmer_repr_t, uint32_t> b) const
		{
			return a.first < b.first;
		}
	} custom_less;

	#pragma omp parallel
	{
	    #pragma omp sections
	    {
	        #pragma omp section
	        {
				kmer_list_father.reserve(total_kmers1);
				read_kmc_db(kmc_db_father, kmer_list_father, kmer_length, kmer_list_mother, total_kmers1 + total_kmers2);
				std::cout << "Finished reading k-mer db (father)." << std::endl;
				std::sort(kmer_list_father.begin(), kmer_list_father.end(), custom_less);
				std::cout << "Finished sorting (father)." << std::endl;
	        }

	        #pragma omp section
	        {
				kmer_list_mother.reserve(total_kmers2);
				read_kmc_db(kmc_db_mother, kmer_list_mother, kmer_length, kmer_list_father, total_kmers1 + total_kmers2);
				std::cout << "Finished reading k-mer db (mother)." << std::endl;
				std::sort(kmer_list_mother.begin(), kmer_list_mother.end(), custom_less);
				std::cout << "Finished sorting (mother)." << std::endl;
	        }
	    }
	}

	uint64_t shared_kmer_count_diff[100] = {0};
	uint64_t shared_kmer_count_diff_low[100] = {0};
	uint64_t intersect_size = 0;
	uint64_t i = 0, j = 0;
	while (true) // parallelize on A..., C..., G..., T...
	{
		if (kmer_list_father[i].first == kmer_list_mother[j].first)
		{
			uint32_t min_val, max_val;
			if (kmer_list_father[i].second < kmer_list_mother[j].second)
			{
				min_val = kmer_list_father[i].second;
				max_val = kmer_list_mother[j].second;
			}
			else
			{
				min_val = kmer_list_mother[j].second;
				max_val = kmer_list_father[i].second;
			}

			uint32_t perc_diff{std::min<uint32_t>((100.0f * max_val / min_val) - 100, 99)};
			if ((i & 0b1111111111111111111111111111) == 0)
				std::cout << min_val << '\t' << max_val << '\t' << perc_diff << std::endl;
			shared_kmer_count_diff[perc_diff]++;
			if (min_val < 5)
				shared_kmer_count_diff_low[perc_diff]++;

			++i, ++j, ++intersect_size;
		}
		else if (kmer_list_father[i].first > kmer_list_mother[j].first)
			++j;
		else
			++i;

		if (i == kmer_list_father.size() || j == kmer_list_mother.size())
			break;
	}

	std::cout << "Father (total/unique/%): " << kmer_list_father.size() << " / "
			  								 << (kmer_list_father.size() - intersect_size) << " / "
											 << (10000 * (kmer_list_father.size() - intersect_size) / kmer_list_father.size()) / 100.0f << "%\n"
	 	  	  << "Mother (total/unique/%): " << kmer_list_mother.size() << " / "
			  								 << (kmer_list_mother.size() - intersect_size) << " / "
 											 << (10000 * (kmer_list_mother.size() - intersect_size) / kmer_list_mother.size()) / 100.0f << "%\n";

	 std::cout << "Stats on shared k-mers (nbr. of k-mers that have more less than X % difference in their count value)\n";
	 for (uint32_t i = 0; i < 100; ++i)
	 	std::cout << "i = " << i << '\t' << shared_kmer_count_diff[i] << "\t(" << shared_kmer_count_diff_low[i] << ")\n";

	return 0;
}
