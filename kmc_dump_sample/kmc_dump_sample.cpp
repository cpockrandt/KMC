#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include <iomanip>
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

struct paired_sequence_file_output {

    using fields_type = fields<field::SEQ, field::ID, field::QUAL>;
    using formats_type = type_list<format_fasta, format_fastq>;
    using stream_char_type = char;

	sequence_file_output<fields_type, formats_type, stream_char_type> pe1;
	sequence_file_output<fields_type, formats_type, stream_char_type> pe2;

	paired_sequence_file_output(std::filesystem::path && path): pe1(path.string() + std::string{".1.fa"}), // TODO: automatic file ending
																pe2(path.string() + std::string{".2.fa"})
	{
		pe1.options.fasta_letters_per_line = 0;
		pe2.options.fasta_letters_per_line = 0;
	}

	void push_back(auto && r1, auto && r2)
	{
		pe1.push_back(r1);
		pe2.push_back(r2);
	}

};

std::string get_perc(auto && count, auto && total)
{
    std::string perc{std::to_string(truncf((count / static_cast<float>(total))*10000)/100)};
    perc.erase(perc.find_last_not_of('0') + 1, std::string::npos);
    return perc;
}

enum class fo_mo_state : uint8_t { NONE = 0, FO_LAST = 1, MO_LAST = 2 };

int _tmain(int argc, char* argv[])
{
	argument_parser myparser{"Haplotype resolution", argc, argv};

	myparser.info.author = "Christopher Pockrandt";
	myparser.info.short_description = "Haplotype resolution of Trios using k-mer counting.";
	myparser.info.version = "0.0.1";

    float perc = 1;
	uint32_t min = 3;
	uint32_t threads = omp_get_max_threads();
	std::filesystem::path kmc_father_path{}, kmc_mother_path{}, reads_child_path1{}, reads_child_path2{}, out_dir{};

	myparser.add_positional_option(kmc_father_path, "Please provide the KMC file of the father.");
	myparser.add_positional_option(kmc_mother_path, "Please provide the KMC file of the mother.");
	myparser.add_positional_option(reads_child_path1, "Please provide the read file of the child (1. PE).");
	myparser.add_positional_option(reads_child_path2, "Please provide the read file of the child (2. PE).");

	myparser.add_option(min, 'm', "min", "Min k-mers (occurring less than X times)", option_spec::DEFAULT);
	myparser.add_option(out_dir, 'o', "out", "Directory to output binned reads.", option_spec::DEFAULT, output_directory_validator{});
	myparser.add_option(threads, 't', "threads", "Number of threads", option_spec::DEFAULT, arithmetic_range_validator{1, 64});
	myparser.add_option(perc, 'p', "perc", "Bin reads where every k-mer in a read has more than p/100% occurrences in a parent.", option_spec::DEFAULT);

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
    	uint32 _mode, _kmer_length2, _counter_size, _lut_prefix_length, _signature_len, _min_count;
    	uint64 _max_count, _total_kmers;
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
	uint64_t global_common = 0;
	uint64_t global_father_only = 0;
	uint64_t global_mother_only = 0;
	uint64_t global_mixed_single_switch = 0;
	uint64_t global_mixed_multiple_switches = 0;

	sequence_file_input fin1{reads_child_path1}, fin2{reads_child_path2};
	auto chunks1 = fin1 // | std::view::transform([] (auto record) { return get<field::SEQ>(record); } )
				           | ranges::view::chunk(500000);
	auto chunks2 = fin2 // | std::view::transform([] (auto record) { return get<field::SEQ>(record); } )
				           | ranges::view::chunk(500000);

	std::vector<typename sequence_file_input<>::record_type> chunks_container1, chunks_container2;
	chunks_container1.reserve(500000);
	chunks_container2.reserve(500000);

	paired_sequence_file_output fout_common{out_dir / "common"},
                                fout_father_only{out_dir / "father_only"},
                                fout_mother_only{out_dir / "mother_only"},
                                fout_mixed_single_switch{out_dir / "mixed_single_switch"},
                                fout_mixed_multiple_switches{out_dir / "mixed_multiple_switches"};
	// fout_father_perc{out_dir / ("father_perc." + std::to_string(perc))}, fout_mother_perc{out_dir / ("mother_perc." + std::to_string(perc))}

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
                	auto chunk = *it1;
                	for (auto chunk_it = chunk.begin(); chunk_it != chunk.end(); chunk_it++)
                	{
                	       chunks_container1.emplace_back(*chunk_it);
                	}
            	}

            	#pragma omp section
            	{
                	chunks_container2.clear();
                	auto chunk = *it2;
                	for (auto chunk_it = chunk.begin(); chunk_it != chunk.end(); chunk_it++)
                	{
                	       chunks_container2.emplace_back(*chunk_it);
                	}
            	}
        	}
    	}

    	std::vector<std::vector<uint32>> counters_father1(threads),
                                         counters_mother1(threads),
                                         counters_father2(threads),
                                         counters_mother2(threads);
        for (uint32_t i = 0; i < threads; ++i)
        {
            counters_father1[i].reserve(250 - kmer_length + 1);
            counters_mother1[i].reserve(250 - kmer_length + 1);
            counters_father2[i].reserve(250 - kmer_length + 1);
            counters_mother2[i].reserve(250 - kmer_length + 1);
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

            auto & counter_father1 = counters_father1[omp_get_thread_num()];
            auto & counter_mother1 = counters_mother1[omp_get_thread_num()];
            auto & counter_father2 = counters_father2[omp_get_thread_num()];
            auto & counter_mother2 = counters_mother2[omp_get_thread_num()];

            counter_father1.clear(); // TODO: necessary?
            counter_mother1.clear(); // TODO: necessary?
            counter_father2.clear(); // TODO: necessary?
            counter_mother2.clear(); // TODO: necessary?

        	kmc_db_father.GetCountersForRead(read1_str, counter_father1);
        	kmc_db_mother.GetCountersForRead(read1_str, counter_mother1);
        	kmc_db_father.GetCountersForRead(read2_str, counter_father2);
        	kmc_db_mother.GetCountersForRead(read2_str, counter_mother2);

        	// std::vector<int32_t> diff;
        	// std::transform(counter_father1.begin(), counter_father1.end(), counter_mother1.begin(), std::back_inserter(diff), std::minus<int>());
        	// debug_stream << counter_father1 << '\n'
        	// 			 << counter_mother1 << '\n'
        	// 			 << diff << std::endl;

        	uint32_t mo = 0, fo = 0, none = 0, mixed = 0; //, mo_perc = 0, fo_perc = 0;
            uint32_t fo_mo_switches = 0;
            fo_mo_state last_fo_mo{fo_mo_state::NONE};
        	for (uint32_t i = 0; i < counter_father1.size(); ++i)
        	{
                // remove errors
                if (counter_father1[i] < 3)
                    counter_father1[i] = 0;
                if (counter_mother1[i] < 3)
                    counter_mother1[i] = 0;

            	if (counter_father1[i] == 0 && counter_mother1[i] >= min)
                {
                    ++mo;
                    if (last_fo_mo == fo_mo_state::FO_LAST)
                    {
                        ++fo_mo_switches;
                    }
                    last_fo_mo = fo_mo_state::MO_LAST;
                }
            	else if (counter_father1[i] >= min && counter_mother1[i] == 0) // TODO: parameter for 0. sufficiently small should be okay!
                {
                    ++fo;
                    if (last_fo_mo == fo_mo_state::MO_LAST)
                    {
                        ++fo_mo_switches;
                    }
                    last_fo_mo = fo_mo_state::FO_LAST;
                }
                else if (counter_father1[i] == 0 && counter_mother1[i] == 0)
                    ++none;
                else
                    ++mixed; // TODO: mixed is not really needed.

            	// if (counter_father1[i] >= perc * counter_mother1[i]) // TODO
            	//    ++fo_perc;
            	// else if (counter_mother1[i] >= perc * counter_father1[i])
            	//    ++mo_perc;
        	}
        	for (uint32_t i = 0; i < counter_father2.size(); ++i)
        	{
                // remove errors
                if (counter_father2[i] < 3)
                    counter_father2[i] = 0;
                if (counter_mother2[i] < 3)
                    counter_mother2[i] = 0;

            	if (counter_father2[i] == 0 && counter_mother2[i] >= min)
                {
            	    ++mo;
                    if (last_fo_mo == fo_mo_state::FO_LAST)
                    {
                        ++fo_mo_switches;
                    }
                    last_fo_mo = fo_mo_state::MO_LAST;
                }
            	else if (counter_father2[i] >= min && counter_mother2[i] == 0)
                {
            	    ++fo;
                    if (last_fo_mo == fo_mo_state::MO_LAST)
                    {
                        ++fo_mo_switches;
                    }
                    last_fo_mo = fo_mo_state::FO_LAST;
                }
               else if (counter_father2[i] == 0 && counter_mother2[i] == 0)
                   ++none;
               else
                   ++mixed;

            	// if (counter_father2[i] >= perc * counter_mother2[i])
            	//    ++fo_perc;
            	// else if (counter_mother2[i] >= perc * counter_father2[i])
            	//    ++mo_perc;
        	}

            if (mo > 0 && fo == 0)
            {
                #pragma omp critical(mother_only)
                {
                    ++global_mother_only;
                	fout_mother_only.push_back(read1, read2);
                }
            }
            else if (mo == 0 && fo > 0)
            {
                #pragma omp critical(father_only)
                {
                    ++global_father_only;
                	fout_father_only.push_back(read1, read2);
                }
            }
            else if (mo == 0 && fo == 0)
            {
                #pragma omp critical(common)
                {
                    ++global_common;
                	fout_common.push_back(read1, read2);
                }
            }
            else // if (mo > 0 && fo > 0)
            {
                if (fo_mo_switches == 1)
                {
                    #pragma omp critical(mixed_single_switch)
                    {
                        ++global_mixed_single_switch;
                    	fout_mixed_single_switch.push_back(read1, read2);

                    	fo_mo[fo][mo]++;
                    }
                }
                else
                {
                    #pragma omp critical(mixed_multiple_switches)
                    {
                        ++global_mixed_multiple_switches;
                    	fout_mixed_multiple_switches.push_back(read1, read2);
                    }
                }
            }

        	#pragma omp atomic
        	++no_reads;

        	// if ((no_reads & 0b11111111111111111) == 0b11111111111111111)
        	// 	std::cout << '.' << std::flush;
    	}
    	std::cout << '.' << std::flush;
	}

    std::cout << std::fixed << std::setprecision(3) << a;
    std::cout << "\n\nTotal:\t" << no_reads << '\n';
    std::cout << "----------------------------------\n";
	std::cout << "Common:\t" << global_common << " (" << (100.0f * global_common / no_reads) << " %)\n";
	std::cout << "Father only:\t" << global_father_only << " (" << (100.0f * global_father_only / no_reads) << " %)\n";
	std::cout << "Mother only:\t" << global_mother_only << " (" << (100.0f * global_mother_only / no_reads) << " %)\n";
	std::cout << "Mixed single switch:\t" << global_mixed_single_switch << " (" << (100.0f * global_mixed_single_switch / no_reads) << " %)\n";
	std::cout << "Mixed multiple switches:\t" << global_mixed_multiple_switches << " (" << (100.0f * global_mixed_multiple_switches / no_reads) << " %)\n";

	// determine last row that contains non-zero values
	uint32_t last_row = 0;
	for (uint32_t j = 0; j < 2 * (250 - kmer_length + 1) + 1; ++j)
	{
    	for (uint32_t i = 0; i < 2 * (250 - kmer_length + 1) + 1; ++i)
    	{
        	if (fo_mo[i][j] > 0 && i > last_row)
                last_row = i;
    	}
	}
	// std::cout << "Last row: " << last_row << '\n';

	for (uint32_t i = 0; i < last_row + 1; ++i)
	{
    	uint32_t last_col = 0;
    	for (uint32_t j = 0; j < 2 * (250 - kmer_length + 1) + 1; ++j)
    	{
        	if (fo_mo[i][j] > 0 && j > last_col)
                last_col = j;
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

	// for (uint32_t i = 1; i < 2 * (250 - kmer_length + 1) + 1; ++i)
	// {
    // 	fa_only += fo_mo[i][0];
    // 	mo_only += fo_mo[0][i];
    //
    // 	for (uint32_t j = 1; j < 2 * (250 - kmer_length + 1) + 1; ++j)
    // 	{
    // 	mixed += fo_mo[i][j];
    // 	}
	// }
	// std::cout << "\n\nTotal:\t" << no_reads << '\n';
	// std::cout << "None:\t" << fo_mo[0][0] << '\n';
	// std::cout << "Father only:\t" << fa_only << '\n';
	// std::cout << "Mother only:\t" << mo_only << '\n';
	// std::cout << "Mixed:\t" << mixed << '\n';

	return 0;
}
