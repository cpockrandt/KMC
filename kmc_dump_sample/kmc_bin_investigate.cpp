#include "omp.h"
#include <iostream>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>

#include "common.hpp"

using namespace seqan3;

template <typename output_t>
inline void run(std::filesystem::path const & p1, std::filesystem::path & p2/*, std::filesystem::path const & out_dir,
                auto & kmc_db, auto & kmc_father, auto & kmc_mother,
                bool const output_details*/, int const threads, float const perc)
{
    std::cout << std::fixed << std::setprecision(2);

    constexpr bool paired_end{false/*std::is_same_v<output_t, paired_sequence_file_output>*/};

    // uint64_t fo_mo[2 * (250 - kmer_length + 1) + 1][2 * (250 - kmer_length + 1) + 1] = {0};
    uint64_t no_reads = 0;
    uint64_t global_father_only = 0;
    uint64_t global_mother_only = 0;
    uint64_t global_mixed_single_switch = 0;
    uint64_t global_mixed_multiple_switches = 0;

    sequence_file_input fin1{p1/*, seqan3::format_fastq{}*/};
    auto chunks1 = fin1 | ranges::view::chunk(10'000);
    std::vector<typename sequence_file_input<>::record_type> chunks_container1;
    chunks_container1.reserve(10'000);
    auto it1{chunks1.begin()};

    if constexpr (!paired_end)
        p2 = p1;
    sequence_file_input fin2{p2/*, seqan3::format_fastq{}*/};
    auto chunks2 = fin2 | ranges::view::chunk(10'000);
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
            // fo_mo_state last_fo_mo{fo_mo_state::NONE};

            auto & read1{chunks_container1[id]};
            auto & read2{chunks_container2[id]};

            // if (output_details)
            // {
                analyze_read_sliding_win<true>(read1, mo, fo, fo_mo_switches, perc);
                if constexpr (paired_end)
                    analyze_read_sliding_win<true>(read2, mo, fo, fo_mo_switches, perc);
            // }
            // else
            // {
            //     analyze_read<false>(read1, mo, fo/*, fo_mo_switches*/);
            //     if constexpr (paired_end)
            //         analyze_read<false>(read2, mo, fo/*, fo_mo_switches*/);
            // }

            // if (!output_details && mo > 0 && fo > 0 && fo_mo_switches > 1)
            //     break;

            if ((mo > 0 && fo == 0)/* || (mo > 2 && fo <= 1)*/)
            {
                #pragma omp critical(mother_only)
                {
                    ++global_mother_only;
                }
            }
            else if ((mo == 0 && fo > 0)/* || (fo > 2 && mo <= 1)*/)
            {
                #pragma omp critical(father_only)
                {
                    ++global_father_only;
                }
            }
            else if (mo > 0 && fo > 0)
            {
                if (fo_mo_switches == 1)
                {
                    #pragma omp critical(mixed_single_switch)
                    {
                        ++global_mixed_single_switch;
                    }
                }
                else
                {
                    #pragma omp critical(mixed_multiple_switches)
                    {
                        ++global_mixed_multiple_switches;
                    }
                }
            }
            else if (mo == 0 && fo == 0)
            {
                #pragma omp critical(common)
                {
                    // ++global_common;
                }
            }
        }

        no_reads += chunks_container1.size();
        if (no_reads > chunks_container1.size())
            std::cout << "\033[7A";

        std::cout << "\nTotal            :\t" << no_reads << '\n';
        std::cout << "-------------------------------------------\n";
        uint64_t const global_common = no_reads - global_father_only - global_mother_only /*- global_mixed_single_switch*/ - global_mixed_multiple_switches;
        std::cout << "Common           :\t" << global_common << " (" << (100.0f * global_common / no_reads) << " %)\n";
        std::cout << "Father only      :\t" << global_father_only << " (" << (100.0f * global_father_only / no_reads) << " %)\n";
        std::cout << "Mother only      :\t" << global_mother_only << " (" << (100.0f * global_mother_only / no_reads) << " %)\n";
        std::cout << "Single switch    :\t" << global_mixed_single_switch << " (" << (100.0f * global_mixed_single_switch / no_reads) << " %)\n";
        std::cout << "Multiple switches:\t" << global_mixed_multiple_switches << " (" << (100.0f * global_mixed_multiple_switches / no_reads) << " %)" << std::flush;

        if constexpr (paired_end)
            ++it2;
    }
}

int main(int argc, char* argv[])
{
    argument_parser myparser{"Haplotype binning", argc, argv};

    myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Haplotype resolution of Trios using k-mer counting.";
    myparser.info.version = "0.0.1";

    bool output_details{false};
    float perc{.8f};
    int threads{omp_get_max_threads()};
    std::filesystem::path reads_child_path{}, reads_child_path1{}, reads_child_path2{}; //, out_dir{}, kmc_path{}, ;

    // std::filesystem::path kmc_father_path, kmc_mother_path{};

    // myparser.add_positional_option(kmc_path, "Please provide the unified KMC file of the parents.");

    // TODO: make sure at least one/two path(s) are passed.
    myparser.add_option(reads_child_path1, 'x', "pe1", "Please provide the read file of the child (1. PE)", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});
    myparser.add_option(reads_child_path2, 'y', "pe2", "Please provide the read file of the child (2. PE)", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});
    myparser.add_option(reads_child_path, 'z', "sp", "Please provide the read file of the child.", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});

    // myparser.add_option(kmc_father_path, 'f', "kmc_father", "Please provide the KMC db of the father.", option_spec::DEFAULT);
    // myparser.add_option(kmc_mother_path, 'm', "kmc_mother", "Please provide the KMC db of the mother.", option_spec::DEFAULT);

    // myparser.add_option(out_dir, 'o', "out", "Directory to output binned reads.", option_spec::DEFAULT, output_directory_validator{});
    myparser.add_option(threads, 't', "threads", "Number of threads", option_spec::DEFAULT, arithmetic_range_validator{1, 64});
    myparser.add_option(perc, 'p', "perc", "A mother/father signal is defined if at least p percent of k-mers in a window of length k contain unique k-mers belonging to monther/father.", option_spec::DEFAULT);
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

    // CKMCFile kmc_db;
    // if (!kmc_db.OpenForRA(kmc_path))
    // {
    //     std::cerr << "Could not open KMC file.\n";
    //     return 1;
    // }
    //
    // CKMCFile kmc_father, kmc_mother;
    // if (!kmc_father.OpenForRA(kmc_father_path) || !kmc_mother.OpenForRA(kmc_mother_path))
    // {
    //     std::cerr << "Could not open KMC files of parents.\n";
    //     return 1;
    // }

    if (paired_end)
        run<paired_sequence_file_output>(reads_child_path1, reads_child_path2/*, out_dir, kmc_db, kmc_father, kmc_mother, output_details*/, threads, perc);
    else
        run<single_sequence_file_output>(reads_child_path, reads_child_path/*, out_dir, kmc_db, kmc_father, kmc_mother, output_details*/, threads, perc);

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
