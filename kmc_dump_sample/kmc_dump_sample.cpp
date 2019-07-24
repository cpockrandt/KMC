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
//     print_info();
//     return EXIT_FAILURE ;
// }
// else
// {
//     uint32 _kmer_length1;
//     uint32 _mode1;
//
//     {
//         uint32 _counter_size;
//         uint32 _lut_prefix_length;
//         uint32 _signature_len;
//         uint32 _min_count;
//         uint64 _max_count;
//         uint64 _total_kmers;
//         kmer_data_base1.Info(_kmer_length1, _mode1, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
//     }
//     std::cout << _kmer_length1 << "-mers\n";
//
//     if(min_count_to_set)
//     if (!(kmer_data_base1.SetMinCount(min_count_to_set)))
//             return EXIT_FAILURE;
//     if(max_count_to_set)
//     if (!(kmer_data_base1.SetMaxCount(max_count_to_set)))
//             return EXIT_FAILURE;
//
//     uint32 c1;
//     CKmerAPI kmer_object1(_kmer_length1);
//     std::string kmer1;
//
//     uint64_t hist[10000] = {0};
//
//     while (kmer_data_base1.ReadNextKmer(kmer_object1, c1))
//     {
//         kmer_object1.to_string(kmer1);
//
//         if (c1 >= 10000)
//             c1 = 9999;
//
//         hist[c1]++;
//     }
//
//     for (uint32_t i = 0; i < 10000; ++i)
//         std::cout << i << '\t' << hist[i] << '\n';
//
//     // fclose(out_file);
//     kmer_data_base1.Close();
// }

struct paired_sequence_file_output {

    using fields_type = fields<field::SEQ, field::ID, field::QUAL>;
    using formats_type = type_list<format_fastq>; // format_fasta,
    using stream_char_type = char;

    sequence_file_output<fields_type, formats_type, stream_char_type> pe1;
    sequence_file_output<fields_type, formats_type, stream_char_type> pe2;

    paired_sequence_file_output(std::filesystem::path && path): pe1(path.string() + std::string{".1.fq"}), // TODO: automatic file ending
                                                                pe2(path.string() + std::string{".2.fq"})
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

// inline void anaylse_read(auto & read, auto & counter, uint16_t & mo, uint16_t & fo, uint16_t & fo_mo_switches, fo_mo_state & last_fo_mo, bool const output_details, auto & kmc_db, uint32 const kmer_length)
// {
//     auto read_view = get<field::SEQ>(read) | view::to_char;
//     auto & qual = get<field::QUAL>(read);
//     std::string read_str(read_view.begin(), read_view.end());
//     kmc_db.GetCountersForRead(read_str, counter);
//
//     qual.resize(read_str.size() - kmer_length + 1);
//
//     for (uint32_t i = 0; i < read_str.size() - kmer_length + 1; ++i)
//     {
//         if (counter[i] == 2)
//         {
//             ++mo;
//             if (last_fo_mo == fo_mo_state::FO_LAST)
//                 ++fo_mo_switches;
//             last_fo_mo = fo_mo_state::MO_LAST;
//             qual[i].assign_char('2');
//         }
//         else if (counter[i] == 1)
//         {
//             ++fo;
//             if (last_fo_mo == fo_mo_state::MO_LAST)
//                 ++fo_mo_switches;
//             last_fo_mo = fo_mo_state::FO_LAST;
//             qual[i].assign_char('1');
//         }
//         else
//         {
//             qual[i].assign_char('.');
//         }
//     }
// }

int _tmain(int argc, char* argv[])
{
    argument_parser myparser{"Haplotype resolution", argc, argv};

    myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Haplotype resolution of Trios using k-mer counting.";
    myparser.info.version = "0.0.1";

    bool output_details{false};
    float perc{1.0f};
    uint32_t min{3};
    int threads{omp_get_max_threads()};
    std::filesystem::path kmc_path{}, reads_child_path1{}, reads_child_path2{}, out_dir{};

    myparser.add_positional_option(kmc_path, "Please provide the unified KMC file of the parents.");
    myparser.add_positional_option(reads_child_path1, "Please provide the read file of the child (1. PE).");
    myparser.add_positional_option(reads_child_path2, "Please provide the read file of the child (2. PE).");

    myparser.add_option(min, 'm', "min", "Min k-mers (occurring less than X times)", option_spec::DEFAULT);
    myparser.add_option(out_dir, 'o', "out", "Directory to output binned reads.", option_spec::DEFAULT, output_directory_validator{});
    myparser.add_option(threads, 't', "threads", "Number of threads", option_spec::DEFAULT, arithmetic_range_validator{1, 64});
    myparser.add_option(perc, 'p', "perc", "Bin reads where every k-mer in a read has more than p/100% occurrences in a parent.", option_spec::DEFAULT);
    myparser.add_option(output_details, 'd', "details", "Output binning details in quality string.", option_spec::DEFAULT);

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
        std::cerr << "Could not open KMC file.\n";
        return 1;
    }

    uint32 kmer_length;
    {
        uint32 _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
        uint64 _max_count, _total_kmers;
        kmc_db.Info(kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
        std::cout << "k-mer length: " << kmer_length << std::endl;
    }

    uint64_t no_reads = 0;
    // uint64_t fo_mo[2 * (250 - kmer_length + 1) + 1][2 * (250 - kmer_length + 1) + 1] = {0};
    uint64_t global_father_only = 0;
    uint64_t global_mother_only = 0;
    uint64_t global_mixed_single_switch = 0;
    uint64_t global_mixed_multiple_switches = 0;

    sequence_file_input fin1{reads_child_path1}, fin2{reads_child_path2};
    auto chunks1 = fin1 // | std::view::transform([] (auto record) { return get<field::SEQ>(record); } )
                           | ranges::view::chunk(1'000'000);
    auto chunks2 = fin2 // | std::view::transform([] (auto record) { return get<field::SEQ>(record); } )
                           | ranges::view::chunk(1'000'000);

    std::vector<typename sequence_file_input<>::record_type> chunks_container1, chunks_container2;
    chunks_container1.reserve(1'000'000);
    chunks_container2.reserve(1'000'000);

    std::cout << std::fixed << std::setprecision(2);

    paired_sequence_file_output // fout_common{out_dir / "common"},
                                fout_father_only{out_dir / "father_only"},
                                fout_mother_only{out_dir / "mother_only"},
                                fout_mixed_single_switch{out_dir / "mixed_single_switch"},
                                fout_mixed_multiple_switches{out_dir / "mixed_multiple_switches"};

    std::vector<std::vector<uint32>> counters(threads);
    for (int i = 0; i < threads; ++i)
        counters[i].reserve(250 - kmer_length + 1);

    auto it1{chunks1.begin()};
    auto it2{chunks2.begin()};
    for (; it1 != chunks1.end(); ++it1, ++it2)
    {
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

        #pragma omp parallel for num_threads(threads)
        for (uint64_t id = 0; id < chunks_container1.size(); ++id)
        {
            uint16_t mo{0}, fo{0}, fo_mo_switches{0};
            fo_mo_state last_fo_mo{fo_mo_state::NONE};

            auto & read1{chunks_container1[id]}, read2{chunks_container2[id]};
            auto & counter = counters[omp_get_thread_num()];

            // anaylse_read(read1, counter, mo, fo, fo_mo_switches, last_fo_mo, output_details, kmc_db, kmer_length);
            // anaylse_read(read2, counter, mo, fo, fo_mo_switches, last_fo_mo, output_details, kmc_db, kmer_length);
            for (auto read : {&read1, &read2})
            {
                auto read_view = get<field::SEQ>(*read) | view::to_char;
                auto & qual = get<field::QUAL>(*read);
                std::string read_str(read_view.begin(), read_view.end());
                kmc_db.GetCountersForRead(read_str, counter);

                qual.resize(read_str.size() - kmer_length + 1);

                for (uint32_t i = 0; i < read_str.size() - kmer_length + 1; ++i)
                {
                    if (counter[i] == 2)
                    {
                        ++mo;
                        if (last_fo_mo == fo_mo_state::FO_LAST)
                            ++fo_mo_switches;
                        last_fo_mo = fo_mo_state::MO_LAST;
                        qual[i].assign_char('2');
                    }
                    else if (counter[i] == 1)
                    {
                        ++fo;
                        if (last_fo_mo == fo_mo_state::MO_LAST)
                            ++fo_mo_switches;
                        last_fo_mo = fo_mo_state::FO_LAST;
                        qual[i].assign_char('1');
                    }
                    else
                    {
                        qual[i].assign_char('.');
                    }
                }
                for (uint32_t i = read_str.size() - kmer_length + 1; i < output_details * read_str.size(); ++i)
                {
                    qual[i].assign_char('.');
                }

                if (!output_details && mo > 0 && fo > 0 && fo_mo_switches > 1)
                    break;
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
            else if (mo > 0 && fo > 0)
            {
                if (fo_mo_switches == 1)
                {
                    #pragma omp critical(mixed_single_switch)
                    {
                        ++global_mixed_single_switch;
                        fout_mixed_single_switch.push_back(read1, read2);
                        // fo_mo[fo][mo]++;
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
            // else if (mo == 0 && fo == 0)
            // {
            //     // #pragma omp critical(common)
            //     // {
            //     //     ++global_common;
            //     //     // fout_common.push_back(read1, read2);
            //     // }
            // }

            // #pragma omp atomic
            // ++no_reads;
        }

        no_reads += chunks_container1.size();
        if (no_reads > chunks_container1.size())
            std::cout << "\033[7A";

        std::cout << "\nTotal            :\t" << no_reads << '\n';
        std::cout << "-------------------------------------------\n";
        uint64_t const global_common = no_reads - global_father_only - global_mother_only - global_mixed_single_switch - global_mixed_multiple_switches;
        std::cout << "Common           :\t" << global_common << " (" << (100.0f * global_common / no_reads) << " %)\n";
        std::cout << "Father only      :\t" << global_father_only << " (" << (100.0f * global_father_only / no_reads) << " %)\n";
        std::cout << "Mother only      :\t" << global_mother_only << " (" << (100.0f * global_mother_only / no_reads) << " %)\n";
        std::cout << "Single switch    :\t" << global_mixed_single_switch << " (" << (100.0f * global_mixed_single_switch / no_reads) << " %)\n";
        std::cout << "Multiple switches:\t" << global_mixed_multiple_switches << " (" << (100.0f * global_mixed_multiple_switches / no_reads) << " %)" << std::flush;
    }
    std::cout << '\n';

    // determine last row that contains non-zero values
    // uint32_t last_row = 0;
    // for (uint32_t j = 0; j < 2 * (250 - kmer_length + 1) + 1; ++j)
    // {
    //     for (uint32_t i = 0; i < 2 * (250 - kmer_length + 1) + 1; ++i)
    //     {
    //         if (fo_mo[i][j] > 0 && i > last_row)
    //             last_row = i;
    //     }
    // }
    //
    // for (uint32_t i = 0; i < last_row + 1; ++i)
    // {
    //     uint32_t last_col = 0;
    //     for (uint32_t j = 0; j < 2 * (250 - kmer_length + 1) + 1; ++j)
    //     {
    //         if (fo_mo[i][j] > 0 && j > last_col)
    //             last_col = j;
    //     }
    //
    //     for (uint32_t j = 0; j < last_col + 1; ++j)
    //     {
    //         if (fo_mo[i][j] > 0)
    //            std::cout << fo_mo[i][j] << '\t';
    //         else
    //            std::cout << '\t';
    //     }
    //     std::cout << '\n';
    // }

    return 0;
}
