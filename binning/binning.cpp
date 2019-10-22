#include "omp.h"
#include "../kmc_dump_sample/stdafx.h"
#include <iostream>
#include <iomanip>
#include "../kmc_api/kmc_file.h"

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/search/algorithm/all.hpp>

#include <seqan3/range/views/to_char.hpp>
#include <range/v3/view/sliding.hpp>

#include "../kmc_dump_sample/common.hpp"
#include "bwt_binning_common.hpp"

using namespace seqan3;

struct program_options
{
    std::filesystem::path reads_child_path{};
    std::filesystem::path out_path{};
    std::string mode{};
    std::vector<std::filesystem::path> kmc_path{};
    std::filesystem::path bwt_path{};

    uint32_t bwt_kmer{25};
    uint8_t bwt_error{0};
    float perc{1.0};
    float perc_low{0.5};
    int threads{omp_get_max_threads()};

    bool details{};
};

inline void analyse(auto & read_str, auto & counter, auto & kmc_db, uint32 const kmer_length, uint32_t & fo, uint32_t & mo)
{
    kmc_db.GetCountersForRead(read_str, counter);

    for (uint32_t i = 0; i < read_str.size() - kmer_length + 1; ++i)
    {
        if (counter[i] == 2)
            ++mo;
        else if (counter[i] == 1)
            ++fo;
    }
}

inline void analyze_single_window(auto const & counter, uint16_t const from, uint16_t const to, uint16_t & window_m, uint16_t & window_f)
{
    window_m = 0;
    window_f = 0;
    for (uint16_t i = from; i < to + 1; ++i)
    {
        if (counter[i] == 2)
            ++window_m;
        else if (counter[i] == 1)
            ++window_f;
    }
}

// inline void algo_window(program_options const & options, std::vector<CKMCFile> & kmc_db)
// {
//     uint32 kmer_length[7] = {0};
//     for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
//     {
//         uint32 _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
//         uint64 _max_count, _total_kmers;
//         kmc_db[db_id].Info(kmer_length[db_id], _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
//     }
//
//     uint64_t no_reads = 0;
//     uint64_t binned_father[7] = {0};
//     uint64_t binned_mother[7] = {0};
//
//     std::vector<std::vector<uint32>> counters(options.threads);
//     for (int i = 0; i < options.threads; ++i)
//     {
//         counters[i].reserve(15'000);
//     }
//
//     sequence_file_input<sequence_file_input_default_traits_dna> fin{options.reads_child_path};
//     auto chunks = fin | ranges::views::chunk(10'000);
//     std::vector<typename decltype(fin)::record_type> chunks_container;
//     chunks_container.reserve(10'000);
//     auto it{chunks.begin()};
//
//     for (; it != chunks.end(); ++it)
//     {
//         loadReads(chunks_container, it);
//
//         #pragma omp parallel for num_threads(options.threads)
//         for (uint64_t id = 0; id < chunks_container.size(); ++id)
//         {
//             auto & read{chunks_container[id]};
//             auto & counter = counters[omp_get_thread_num()];
//             auto read_view = get<field::SEQ>(read) | views::to_char;
//             std::string read_str(read_view.begin(), read_view.end());
//
//             for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
//             {
//             }
//         }
//
//         no_reads += chunks_container.size();
//         if (no_reads > chunks_container.size())
//            std::cout << "\033[" << (2 + 2 * kmc_db.size() ) << "A";
//
//         std::cout << "Total             :\t" << no_reads << '\n';
//         uint64_t unbinned = no_reads;
//         for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
//         {
//             unbinned -= binned_father[db_id];
//             unbinned -= binned_mother[db_id];
//             std::cout << "Father (" << kmer_length[db_id] << ")       :\t" << binned_father[db_id] << " (" << (100.0f * binned_father[db_id] / no_reads) << " %)\n";
//             std::cout << "Mother (" << kmer_length[db_id] << ")       :\t" << binned_mother[db_id] << " (" << (100.0f * binned_mother[db_id] / no_reads) << " %)\n";
//         }
//         std::cout << "Unbinned          :\t" << unbinned << " (" << (100.0f * unbinned / no_reads) << " %)\n";
//         std::cout.flush();
//     }
// }

inline void write_qualities(auto & qual, auto & counter, uint32 const kmer_length)
{
    for (uint32_t i = 0; i < qual.size() - kmer_length + 1; ++i)
    {
        if (counter[i] == 2)
            qual[i].assign_char('2');
        else if (counter[i] == 1)
            qual[i].assign_char('1');
        else
            qual[i].assign_char('.');
    }
    for (uint32_t i = qual.size() - kmer_length + 1; i < qual.size(); ++i)
        qual[i].assign_char('.');
}

template <bool algo_majority>
inline void algo_kmer(program_options const & options, std::vector<CKMCFile> & kmc_db)
{
    uint32 kmer_length[7] = {0};
    for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
    {
        uint32 _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
        uint64 _max_count, _total_kmers;
        kmc_db[db_id].Info(kmer_length[db_id], _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
    }

    uint64_t no_reads = 0;
    uint64_t binned_father[7] = {0};
    uint64_t binned_mother[7] = {0};

    std::string format{"fa"};
    if (options.details)
        format = "fq";

    sequence_file_output fout_father   {options.out_path / ("father." + format)},
                         fout_mother   {options.out_path / ("mother." + format)},
                         fout_unbinned {options.out_path / ("unbinned." + format)};

    fout_father.options.fasta_letters_per_line   = 0;
    fout_mother.options.fasta_letters_per_line   = 0;
    fout_unbinned.options.fasta_letters_per_line = 0;

    std::vector<std::vector<uint32>> counters(options.threads);
    for (int i = 0; i < options.threads; ++i)
    {
        counters[i].reserve(15'000);
    }

    sequence_file_input<sequence_file_input_default_traits_dna> fin{options.reads_child_path};
    auto chunks = fin | ranges::views::chunk(10'000);
    std::vector<typename decltype(fin)::record_type> chunks_container;
    chunks_container.reserve(10'000);
    auto it{chunks.begin()};

    for (; it != chunks.end(); ++it)
    {
        loadReads(chunks_container, it);

        #pragma omp parallel for num_threads(options.threads)
        for (uint64_t id = 0; id < chunks_container.size(); ++id)
        {
            auto & read{chunks_container[id]};
            auto & counter = counters[omp_get_thread_num()];

            auto read_view = get<field::SEQ>(read) | views::to_char;
            std::string read_str(read_view.begin(), read_view.end());

            if constexpr (algo_majority)
            {
                std::cerr << "Majority mode does not have --detail implemented yet. If you don't care, remove this error message and exit();\n";
                exit(1);
                for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
                {
                    uint32_t fo = 0, mo = 0;
                    analyse(read_str, counter, kmc_db[db_id], kmer_length[db_id], fo, mo);

                    if (fo > mo) // if (fo >= 15 && mo <= 5)
                    {
                        #pragma omp critical(fo_inc)
                        {
                            ++binned_father[db_id];
                            fout_father.push_back(read);
                        }
                        break;
                    }
                    else if (mo > fo) // else if (mo >= 15 && fo <= 5)
                    {
                        #pragma omp critical(mo_inc)
                        {
                            ++binned_mother[db_id];
                            fout_mother.push_back(read);
                        }
                        break;
                    }
                }
            }
            else // algo_window
            {
                std::vector<std::vector<uint32>> counters_for_unbinned;

                bool isBinned = false;
                for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
                {
                    uint32_t mo = 0;
                    uint32_t fo = 0;
                    // fo_mo_switches = 0;
                    uint32_t mo_ambig = 0;
                    uint32_t fo_ambig = 0;
                    // fo_mo_state last_fo_mo{fo_mo_state::NONE};

                    kmc_db[db_id].GetCountersForRead(read_str, counter);

                    uint16_t window_m, window_f;
                    analyze_single_window(counter, 0, kmer_length[db_id] - 1, window_m, window_f);
                    // debug_stream << "init: " << window_f << " ... " << window_m << " ... " << std::endl;

                    for (uint32_t i = 0; i < read_str.size() - kmer_length[db_id] + 2; ++i)
                    {
                        if (window_m >= options.perc * kmer_length[db_id] || window_f >= options.perc * kmer_length[db_id])
                        {
                            if (window_m > window_f)
                            {
                                ++mo;
                                // if (last_fo_mo == fo_mo_state::FO_LAST)
                                //     ++fo_mo_switches;
                                // last_fo_mo = fo_mo_state::MO_LAST;
                            }
                            else
                            {
                                ++fo;
                                // if (last_fo_mo == fo_mo_state::MO_LAST)
                                //     ++fo_mo_switches;
                                // last_fo_mo = fo_mo_state::FO_LAST;
                            }
                            i += kmer_length[db_id] - 1; // i will be incremented by for loop!
                            analyze_single_window(counter, i + 1, i + 1 + kmer_length[db_id] - 1, window_m, window_f);
                            // debug_stream << "jump: " << window_f << " ... " << window_m << " ... " << std::endl;
                            continue;
                        }
                        else if (window_m >= options.perc_low * kmer_length[db_id] || window_f >= options.perc_low * kmer_length[db_id])
                        {
                            if (window_m > window_f)
                                ++mo_ambig;
                            else
                                ++fo_ambig;
                        }

                        if (i < read_str.size() - kmer_length[db_id] + 1)
                        {
                            // sliding window: increment if it is moving into the window
                            if (counter[i + kmer_length[db_id]] == 2)
                                ++window_m;
                            else if (counter[i + kmer_length[db_id]] == 1)
                                ++window_f;

                            // sliding window: decrement if it is moving out of the window
                            if (counter[i] == 2)
                                --window_m;
                            else if (counter[i] == 1)
                                --window_f;
                            // debug_stream << "slide (" << i << "): " << window_f << " ... " << window_m << " ... " << "( " << counter[i + kmer_length[db_id]] << " " << counter[i] << " )" << std::endl;
                        }
                    }

                    std::vector<phred42> qual{};

                    if (fo > 0 && mo == 0 && mo_ambig == 0)
                    {
                        isBinned = true;
                        if (options.details)
                        {
                            qual.resize(read_str.size());
                            write_qualities(qual, counter, kmer_length[db_id]);
                        }

                        #pragma omp critical(fo_inc)
                        {
                            ++binned_father[db_id];
                            if (options.details)
                            {
                                std::string read_id = get<field::ID>(read) + "_k" + std::to_string(kmer_length[db_id]);
                                fout_father.push_back(std::tie(get<field::SEQ>(read), read_id, qual));
                            }
                            else
                                fout_father.push_back(read);
                        }
                        break;
                    }
                    else if (mo > 0 && fo == 0 && fo_ambig == 0)
                    {
                        isBinned = true;
                        if (options.details)
                        {
                            qual.resize(read_str.size());
                            write_qualities(qual, counter, kmer_length[db_id]);
                        }

                        #pragma omp critical(mo_inc)
                        {
                            ++binned_mother[db_id];
                            if (options.details)
                            {
                                std::string read_id = get<field::ID>(read) + "_k" + std::to_string(kmer_length[db_id]);
                                fout_mother.push_back(std::tie(get<field::SEQ>(read), read_id, qual));
                            }
                            else
                                fout_mother.push_back(read);
                        }
                        break;
                    }
                    else
                    {
                        if (options.details)
                            counters_for_unbinned.push_back(counter);
                    }

                    // debug_stream << counter << std::endl << fo << " ... " << mo << " ... " << fo_mo_switches << std::endl;
                    // exit(0);
                }

                if (!isBinned)
                {
                    #pragma omp critical(output_unbinned)
                    {
                        if (options.details)
                        {
                            std::vector<phred42> qual{};
                            qual.resize(read_str.size());

                            for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
                            {
                                std::string read_id = get<field::ID>(read) + "_k" + std::to_string(kmer_length[db_id]);
                                write_qualities(qual, counters_for_unbinned[db_id], kmer_length[db_id]);
                                fout_unbinned.push_back(std::tie(get<field::SEQ>(read), read_id, qual));
                            }
                        }
                        else
                        {
                            fout_unbinned.push_back(read);
                        }
                    }
                }
            }
        }

        no_reads += chunks_container.size();
        if (no_reads > chunks_container.size())
           std::cout << "\033[" << (3 + 2 * kmc_db.size()) << "A";

        std::cout << "Total      :\t" << no_reads << '\n';
        uint64_t unbinned = no_reads;
        for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
        {
            unbinned -= binned_father[db_id];
            unbinned -= binned_mother[db_id];
            std::cout << "Father (" << kmer_length[db_id] << "):\t" << binned_father[db_id] << " (" << (100.0f * binned_father[db_id] / no_reads) << " %)\n";
            std::cout << "Mother (" << kmer_length[db_id] << "):\t" << binned_mother[db_id] << " (" << (100.0f * binned_mother[db_id] / no_reads) << " %)\n";
        }
        std::cout << "Unbinned   :\t" << unbinned << " (" << (100.0f * unbinned / no_reads) << " %)\n";
        std::cout << "Binned     :\t" << (no_reads-unbinned) << " (" << (100.0f * (no_reads-unbinned) / no_reads) << " %)\n";
        std::cout.flush();
    }
}

// inline void algo_bwt(program_options const & options, bwt_binner const & bwt)
// {
//     // configuration const cfg = search_cfg::max_error{search_cfg::total{1},
//     //                                                 search_cfg::substitution{1},
//     //                                                 search_cfg::insertion{0},
//     //                                                 search_cfg::deletion{0}};
//     detail::search_param errors{options.bwt_error,options.bwt_error,0,0}; // {1,1,0,0};
//
//     uint64_t no_reads = 0;
//     uint64_t binned_father = 0;
//     uint64_t binned_mother = 0;
//
//     sequence_file_input<sequence_file_input_default_traits_dna> fin{options.reads_child_path};
//     auto chunks = fin | ranges::views::chunk(100);
//     std::vector<typename decltype(fin)::record_type> chunks_container;
//     chunks_container.reserve(100);
//     auto it{chunks.begin()};
//
//     sdsl::rank_support_il<1, 512> indicator_rank(&bwt.indicator);
//
//     for (; it != chunks.end(); ++it)
//     {
//         loadReads(chunks_container, it);
//
//         #pragma omp parallel for num_threads(options.threads)
//         for (uint64_t id = 0; id < chunks_container.size(); ++id)
//         {
//             auto & read{chunks_container[id]};
//             // auto & counter = counters[omp_get_thread_num()];
//             auto read_view = get<field::SEQ>(read) | std::views::transform([](dna5 c){ return static_cast<dna4>(c); });;
//             // std::string read_str(read_view.begin(), read_view.end());
//
//             uint32_t mo = 0;
//             uint32_t fo = 0;
//
//             auto kmers = read_view | ranges::views::sliding(options.bwt_kmer);
//
//             for (auto const & kmer : kmers)
//             {
//                 std::vector<bi_fm_index_cursor<bi_fm_index<dna4, text_layout::collection, my_sdsl_wt_index_type>>> hits_cursors;
//                 auto delegate = [&hits_cursors](auto const & cursor)
//                 {
//                     hits_cursors.push_back(cursor);
//                 };
//
//                 // auto hits = search(kmer, bwt.index, cfg) << '\n';
//                 detail::search_algo_bi<false>(bwt.index, kmer, errors, delegate);
//
//                 // uint64_t total_hits = 0;
//                 uint64_t mother_hits = 0;
//                 uint64_t father_hits = 0;
//                 for (auto const & cursor : hits_cursors)
//                 {
//                     mother_hits += indicator_rank(cursor.fwd_rb + 1) - indicator_rank(cursor.fwd_lb);
//                     // total_hits += 1 + cursor.fwd_rb - cursor.fwd_lb;
//                     father_hits += (1 + cursor.fwd_rb - cursor.fwd_lb) - mother_hits;
//                 }
//                 if (father_hits > 0 && mother_hits == 0)
//                     ++fo;
//                 else if (mother_hits > 0 && father_hits == 0)
//                     ++mo;
//
//                 // debug_stream << hits_cursors << std::endl;
//             }
//
//             if (fo >= 15 && mo <= 5)
//             {
//               #pragma omp critical(fo_inc)
//                 ++binned_father;
//             }
//             else if (mo >= 15 && fo <= 5)
//             {
//               #pragma omp critical(mo_inc)
//                 ++binned_mother;
//             }
//         }
//
//         no_reads += chunks_container.size();
//         if (no_reads > chunks_container.size())
//            std::cout << "\033[" << (3 + 2 * 1 ) << "A";
//
//         std::cout << "Total     :\t" << no_reads << '\n';
//         uint64_t unbinned = no_reads;
//         unbinned -= binned_father;
//         unbinned -= binned_mother;
//         std::cout << "Father    :\t" << binned_father << " (" << (100.0f * binned_father / no_reads) << " %)\n";
//         std::cout << "Mother    :\t" << binned_mother << " (" << (100.0f * binned_mother / no_reads) << " %)\n";
//         std::cout << "Unbinned  :\t" << unbinned << " (" << (100.0f * unbinned / no_reads) << " %)\n";
//         std::cout << "Binned    :\t" << (no_reads-unbinned) << " (" << (100.0f * (no_reads-unbinned) / no_reads) << " %)\n";
//         std::cout.flush();
//     }
// }

inline int run(program_options const & options)
{
    if (options.mode == "window" || options.mode == "majority")
    {
        std::vector<CKMCFile> kmc_db;
        kmc_db.resize(options.kmc_path.size());
        for (uint32_t i = 0; i < options.kmc_path.size(); ++i)
        {
            if (!kmc_db[i].OpenForRA(options.kmc_path[i]))
            {
                std::cerr << "Could not open KMC file.\n";
                return 1;
            }
        }

        if (options.mode == "majority")
        {
            algo_kmer<true>(options, kmc_db);
        }
        else if (options.mode == "window")
        {
            algo_kmer<false>(options, kmc_db);
        }
    }
    // else if (options.mode == "bwt")
    // {
    //     bwt_binner bwt;
    //     {
    //         std::ifstream is{options.bwt_path, std::ios::binary};
    //         cereal::BinaryInputArchive iarchive{is};
    //         iarchive(bwt);
    //     }
    //     algo_bwt(options, bwt);
    // }

    return 0;

    // std::cout << "\n\n";
    // determine last row that contains non-zero values
    // uint32_t last_row = 0;
    // for (uint32_t j = 0; j < 2 * 250 + 1; ++j)
    // {
    //     for (uint32_t i = 0; i < 2 * 250 + 1; ++i)
    //     {
    //         if (fo_mo[i][j] > 0 && i > last_row)
    //             last_row = i;
    //     }
    // }
    //
    // for (uint32_t i = 0; i < last_row + 1; ++i)
    // {
    //     uint32_t last_col = 0;
    //     for (uint32_t j = 0; j < 2 * 250 + 1; ++j)
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

int _tmain(int argc, char* argv[])
{
    argument_parser parser{"Haplotype-binning", argc, argv};

    parser.info.author = "Christopher Pockrandt";
    parser.info.short_description = "Haplotype resolution of Trios using k-mer counting.";
    parser.info.version = "0.0.1";

    program_options options;

    parser.add_option(options.reads_child_path, 'z', "sp", "Please provide the read file of the child.", option_spec::REQUIRED, input_file_validator{{"fa", "fasta", "fq","fastq","gz"}});

    parser.add_option(options.out_path, 'o', "out", "Please provide the path to the output read files.", option_spec::DEFAULT);

    // NOTE: when using multiple k-mer databases it outputs the k-mers for the K with which it was binned.
    parser.add_flag(options.details, 'd', "details", "Output read files as fastq with qualities indicating unique k-mers.");

    parser.add_option(options.mode, 'm', "mode", "Select mode to run.", option_spec::REQUIRED, value_list_validator{{"window", "majority", "bwt"}});

    parser.add_option(options.kmc_path, 'k', "kmc-db", "Unified KMC file of parents (can give multiple, shorter k-mers first).");
    // parser.add_option(options.bwt_path, 'b', "bwt", "Path to BWT index.");

    parser.add_option(options.perc, 'p', "perc", "Percentage of k-mers in sliding window (only for window mode)", option_spec::DEFAULT, arithmetic_range_validator{0, 1});
    parser.add_option(options.perc_low, 'q', "perc-lower", "Percentage of k-mers in sliding window to detect ambiguouity (only for window mode)", option_spec::DEFAULT, arithmetic_range_validator{0, 1});
    // parser.add_option(options.bwt_kmer, 'l', "bwt-kmer", "k-mer size for BWT", option_spec::DEFAULT, arithmetic_range_validator{10, 250});
    // parser.add_option(options.bwt_error, 'e', "bwt-error", "k-mer size for BWT", option_spec::DEFAULT, arithmetic_range_validator{0, 2});

    parser.add_option(options.threads, 't', "threads", "Number of threads", option_spec::DEFAULT, arithmetic_range_validator{1, 64});

    try
    {
        parser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    std::cout << '\n' << std::fixed << std::setprecision(2);

    return run(options);
}
