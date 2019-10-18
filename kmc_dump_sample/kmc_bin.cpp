#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include "../kmc_api/kmc_file.h"

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>

#include "common.hpp"

using namespace seqan3;

template <bool output_details>
inline void anaylse_read(auto & read, auto & counter, uint16_t & mo, uint16_t & fo, uint16_t & fo_mo_switches, fo_mo_state & last_fo_mo, auto & kmc_db, uint32 const kmer_length)
{
    auto read_view = get<field::SEQ>(read) | views::to_char;
    auto & qual = get<field::QUAL>(read);
    std::string read_str(read_view.begin(), read_view.end());
    kmc_db.GetCountersForRead(read_str, counter);

    if (output_details)
        qual.resize(read_str.size());

    for (uint32_t i = 0; i < read_str.size() - kmer_length + 1; ++i)
    {
        if (counter[i] == 2)
        {
            ++mo;
            if (last_fo_mo == fo_mo_state::FO_LAST)
                ++fo_mo_switches;
            last_fo_mo = fo_mo_state::MO_LAST;
            if (output_details)
                qual[i].assign_char('2');
        }
        else if (counter[i] == 1)
        {
            ++fo;
            if (last_fo_mo == fo_mo_state::MO_LAST)
                ++fo_mo_switches;
            last_fo_mo = fo_mo_state::FO_LAST;
            if (output_details)
                qual[i].assign_char('1');
        }
        else if (output_details)
        {
            qual[i].assign_char('.');
        }
    }

    for (uint32_t i = read_str.size() - kmer_length + 1; i < output_details * read_str.size(); ++i)
        qual[i].assign_char('.');
}

template <typename output_t>
inline void run(std::filesystem::path const & p1, std::filesystem::path & p2, std::filesystem::path const & out_dir,
                auto & kmc_db, auto & kmc_father, auto & kmc_mother,
                bool const output_details, int const threads)
{
    uint32 kmer_length;
    {
        uint32 _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
        uint64 _max_count, _total_kmers;
        kmc_db.Info(kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
        std::cout << "k-mer length: " << kmer_length << std::endl;
    }

    std::cout << std::fixed << std::setprecision(2);

    constexpr bool paired_end{std::is_same_v<output_t, paired_sequence_file_output>};

    // uint64_t fo_mo[2 * (250 - kmer_length + 1) + 1][2 * (250 - kmer_length + 1) + 1] = {0};
    uint64_t no_reads = 0;
    uint64_t global_father_only = 0;
    uint64_t global_mother_only = 0;
    uint64_t global_mixed_single_switch = 0;
    uint64_t global_mixed_multiple_switches = 0;

    output_t fout_common{out_dir / "common"},
             fout_father_only{out_dir / "father_only"},
             fout_mother_only{out_dir / "mother_only"},
             fout_mixed_single_switch{out_dir / "mixed_single_switch"},
             fout_mixed_multiple_switches{out_dir / "mixed_multiple_switches"};

    std::vector<std::vector<uint32>> counters(threads);
    std::vector<std::vector<uint32>> counters1(threads), counters2(threads);
    for (int i = 0; i < threads; ++i)
    {
        counters[i].reserve(250 - kmer_length + 1);
        counters1[i].reserve(250 - kmer_length + 1);
        counters2[i].reserve(250 - kmer_length + 1);
    }

    sequence_file_input fin1{p1};
    auto chunks1 = fin1 | ranges::views::chunk(10'000);
    std::vector<typename sequence_file_input<>::record_type> chunks_container1;
    chunks_container1.reserve(10'000);
    auto it1{chunks1.begin()};

    if constexpr (!paired_end)
        p2 = p1;
    sequence_file_input fin2{p2};
    auto chunks2 = fin2 | ranges::views::chunk(10'000);
    std::vector<typename sequence_file_input<>::record_type> chunks_container2;
    chunks_container2.reserve(10'000);
    auto it2{chunks2.begin()};

    for (; it1 != chunks1.end(); ++it1)
    {
        if (paired_end)
            loadReads(chunks_container1, chunks_container2, it1, it2);
        else
            loadReads(chunks_container1, it1);

        #pragma omp parallel for num_threads(threads)
        for (uint64_t id = 0; id < chunks_container1.size(); ++id)
        {
            uint16_t mo{0}, fo{0}, fo_mo_switches{0};
            fo_mo_state last_fo_mo{fo_mo_state::NONE};

            auto & read1{chunks_container1[id]};
            auto & read2{chunks_container2[id]};
            auto & counter = counters[omp_get_thread_num()];

            if (output_details)
            {
                anaylse_read<true>(read1, counter, mo, fo, fo_mo_switches, last_fo_mo, kmc_db, kmer_length);
                if constexpr (paired_end)
                    anaylse_read<true>(read2, counter, mo, fo, fo_mo_switches, last_fo_mo, kmc_db, kmer_length);
            }
            else
            {
                anaylse_read<false>(read1, counter, mo, fo, fo_mo_switches, last_fo_mo, kmc_db, kmer_length);
                if constexpr (paired_end)
                    anaylse_read<false>(read2, counter, mo, fo, fo_mo_switches, last_fo_mo, kmc_db, kmer_length);
            }

            // if (!output_details && mo > 0 && fo > 0 && fo_mo_switches > 1)
            //     break;

            if ((mo > 0 && fo == 0)/* || (mo > 2 && fo <= 1)*/)
            {
                #pragma omp critical(mother_only)
                {
                    ++global_mother_only;
                    // fout_mother_only.push_back(read1, read2);
                }
            }
            else if ((mo == 0 && fo > 0)/* || (fo > 2 && mo <= 1)*/)
            {
                #pragma omp critical(father_only)
                {
                    ++global_father_only;
                    // fout_father_only.push_back(read1, read2);
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
            // else
            // {
            //     auto & counter1 = counters1[omp_get_thread_num()];
            //     auto & counter2 = counters2[omp_get_thread_num()];
            //
            //     auto read_view = get<field::SEQ>(read1) | views::to_char;
            //     std::string read_str(read_view.begin(), read_view.end());
            //
            //     kmc_father.GetCountersForRead(read_str, counter1);
            //     kmc_mother.GetCountersForRead(read_str, counter2);
            //
            //     for (uint64_t x = 0; x < counter1.size(); ++x)
            //     {
            //         if (counter1[x] < 10 || counter2[x] < 10)
            //             if (counter1[x] != 0 && counter2[x] != 0)
            //                 debug_stream << counter1[x] << '\t' << counter2[x] << '\n';
            //     }
            //     debug_stream << "---------------------------------------\n";
            // }
            else if (mo == 0 && fo == 0)
            {
                #pragma omp critical(common)
                {
                    // ++global_common;
                    // fout_common.push_back(read1, read2);
                }
            }
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

        if constexpr (paired_end)
            ++it2;
    }
}

int _tmain(int argc, char* argv[])
{
    argument_parser myparser{"Haplotype binning", argc, argv};

    myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Haplotype resolution of Trios using k-mer counting.";
    myparser.info.version = "0.0.1";

    bool output_details{false};
    float perc{1.0f};
    int threads{omp_get_max_threads()};
    std::filesystem::path kmc_path{}, reads_child_path{}, reads_child_path1{}, reads_child_path2{}, out_dir{};

    std::filesystem::path kmc_father_path, kmc_mother_path{};

    myparser.add_positional_option(kmc_path, "Please provide the unified KMC file of the parents.");

    // TODO: make sure at least one/two path(s) are passed.
    myparser.add_option(reads_child_path1, 'x', "pe1", "Please provide the read file of the child (1. PE)", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});
    myparser.add_option(reads_child_path2, 'y', "pe2", "Please provide the read file of the child (2. PE)", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});
    myparser.add_option(reads_child_path, 'z', "sp", "Please provide the read file of the child.", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});

    myparser.add_option(kmc_father_path, 'f', "kmc_father", "Please provide the KMC db of the father.", option_spec::DEFAULT);
    myparser.add_option(kmc_mother_path, 'm', "kmc_mother", "Please provide the KMC db of the mother.", option_spec::DEFAULT);

    myparser.add_option(out_dir, 'o', "out", "Directory to output binned reads.", option_spec::DEFAULT, output_directory_validator{});
    myparser.add_option(threads, 't', "threads", "Number of threads", option_spec::DEFAULT, arithmetic_range_validator{1, 64});
    myparser.add_option(perc, 'p', "perc", "Bin reads where every k-mer in a read has more than p/100% occurrences in a parent.", option_spec::DEFAULT);
    myparser.add_flag(output_details, 'd', "details", "Output binning details in quality string.", option_spec::DEFAULT);

    try
    {
        myparser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    if (reads_child_path1.empty() != reads_child_path2.empty())
    {
        std::cerr << "Cannot specify only one paired-end file. Use --sp for single-end files instead.\n";
        return 1;
    }

    if (reads_child_path.empty() && reads_child_path1.empty())
    {
        std::cerr << "No read file specified.\n";
        return 1;
    }

    if (!reads_child_path.empty() && !reads_child_path1.empty())
    {
        std::cerr << "You cannot specify both, paired-end and single-end read files.\n";
        return 1;
    }

    bool const paired_end{reads_child_path.empty()};

    CKMCFile kmc_db;
    if (!kmc_db.OpenForRA(kmc_path))
    {
        std::cerr << "Could not open KMC file.\n";
        return 1;
    }

    CKMCFile kmc_father, kmc_mother;
    if (!kmc_father.OpenForRA(kmc_father_path) || !kmc_mother.OpenForRA(kmc_mother_path))
    {
        std::cerr << "Could not open KMC files of parents.\n";
        return 1;
    }

    if (paired_end)
        run<paired_sequence_file_output>(reads_child_path1, reads_child_path2, out_dir, kmc_db, kmc_father, kmc_mother, output_details, threads);
    else
        run<single_sequence_file_output>(reads_child_path, reads_child_path, out_dir, kmc_db, kmc_father, kmc_mother, output_details, threads);

    std::cout << '\n';
    return 0;

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
}
