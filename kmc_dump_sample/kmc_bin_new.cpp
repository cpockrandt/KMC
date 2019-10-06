#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include "../kmc_api/kmc_file.h"

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>

#include <seqan3/range/views/to_char.hpp>

#include "common.hpp"

using namespace seqan3;

inline void analyse(auto & read, auto & counter, auto & kmc_db, uint32 const kmer_length, uint32_t & fo, uint32_t & mo)
{
    auto read_view = get<field::SEQ>(read) | views::to_char;
    std::string read_str(read_view.begin(), read_view.end());
    kmc_db.GetCountersForRead(read_str, counter);

    for (uint32_t i = 0; i < read_str.size() - kmer_length + 1; ++i)
    {
        if (counter[i] == 2)
            ++mo;
        else if (counter[i] == 1)
            ++fo;
    }
}

template <typename output_t>
inline void run(std::filesystem::path const & p1, std::filesystem::path & p2,
                auto & kmc_db,/* auto & kmc_father, auto & kmc_mother,*/ std::filesystem::path const & out,
                int const threads)
{
    uint32 kmer_length[7] = {0};
    for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
    {
        uint32 _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
        uint64 _max_count, _total_kmers;
        kmc_db[db_id].Info(kmer_length[db_id], _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
        // std::cout << "k-mer length: " << kmer_length[db_id] << std::endl;
    }

    std::cout << '\n' << std::fixed << std::setprecision(2);

    constexpr bool paired_end{std::is_same_v<output_t, paired_sequence_file_output>};
    // std::string p_out{p1};
    // p_out.insert(static_cast<std::string>(p1).find_last_of('.'), ".qual" + std::to_string(kmer_length));
    //std::string out_str{out};
    //output_t fout_all{out_str.substr(0, out_str.find_last_of('.'))};

    uint64_t no_reads = 0;
    uint64_t binned_father[7] = {0};
    uint64_t binned_mother[7] = {0};
    uint64_t fo_mo[2 * 250 + 1][2 * 250 + 1] = {0};
    // uint64_t fo_mo[2 * (250 - kmer_length + 1) + 1][2 * (250 - kmer_length + 1) + 1] = {0};

    std::vector<std::vector<uint32>> counters(threads);
    for (int i = 0; i < threads; ++i)
    {
        counters[i].reserve(15'000);
    }

    sequence_file_input<sequence_file_input_default_traits_dna> fin1{p1};
    auto chunks1 = fin1 | ranges::view::chunk(10'000);
    std::vector<typename decltype(fin1)::record_type> chunks_container1;
    chunks_container1.reserve(10'000);
    auto it1{chunks1.begin()};

    if constexpr (!paired_end)
        p2 = p1;
    sequence_file_input<sequence_file_input_default_traits_dna> fin2{p2};
    auto chunks2 = fin2 | ranges::view::chunk(10'000);
    std::vector<typename decltype(fin2)::record_type> chunks_container2;
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
            auto & read1{chunks_container1[id]};
            auto & read2{chunks_container2[id]};
            auto & counter = counters[omp_get_thread_num()];

            for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
            {
                uint32_t fo = 0, mo = 0;
                analyse(read1, counter, kmc_db[db_id], kmer_length[db_id], fo, mo);
                if constexpr (paired_end)
                    analyse(read2, counter, kmc_db[db_id], kmer_length[db_id], fo, mo);

                // if (fo >= 15 && mo <= 5)
                // {
            	//     #pragma omp critical(fo_inc)
                //     ++binned_father[db_id];
                //     break;
                // }
                // else if (mo >= 15 && fo <= 5)
                // {
            	//     #pragma omp critical(mo_inc)
                //     ++binned_mother[db_id];
                //     break;
                // }

        	    #pragma omp critical
        	    fo_mo[std::min<uint32_t>(fo, 500)][std::min<uint32_t>(mo, 500)]++;
            }

            // #pragma omp critical(common)
            // {
            //     fout_all.push_back(read1, read2);
            // }
        }

        no_reads += chunks_container1.size();
        if (no_reads > chunks_container1.size())
           std::cout << "\033[" << (2 + 2*kmc_db.size() ) << "A";

        std::cout << "Total            :\t" << no_reads << '\n';
        // uint64_t unbinned = no_reads;
        // for (uint32_t db_id = 0; db_id < kmc_db.size(); ++db_id)
        // {
        //     unbinned -= binned_father[db_id];
        //     unbinned -= binned_mother[db_id];
        //     std::cout << "Father (" << kmer_length[db_id] << ")       :\t" << binned_father[db_id] << " (" << (100.0f * binned_father[db_id] / no_reads) << " %)\n";
        //     std::cout << "Mother (" << kmer_length[db_id] << ")       :\t" << binned_mother[db_id] << " (" << (100.0f * binned_mother[db_id] / no_reads) << " %)\n";
        // }
        // std::cout << "Unbinned         :\t" << unbinned << " (" << (100.0f * unbinned / no_reads) << " %)\n";
        std::cout.flush();

        if constexpr (paired_end)
            ++it2;
    }

std::cout << "\n\n";

     // determine last row that contains non-zero values
     uint32_t last_row = 0;
     for (uint32_t j = 0; j < 2 * 250 + 1; ++j)
     {
         for (uint32_t i = 0; i < 2 * 250 + 1; ++i)
         {
             if (fo_mo[i][j] > 0 && i > last_row)
                 last_row = i;
         }
     }

     for (uint32_t i = 0; i < last_row + 1; ++i)
     {
         uint32_t last_col = 0;
         for (uint32_t j = 0; j < 2 * 250 + 1; ++j)
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
}

int _tmain(int argc, char* argv[])
{
    argument_parser myparser{"Haplotype-binning", argc, argv};

    myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Haplotype resolution of Trios using k-mer counting.";
    myparser.info.version = "0.0.1";

    int threads{omp_get_max_threads()};
    std::filesystem::path reads_child_path{}, reads_child_path1{}, reads_child_path2{}, out_path{};
    std::vector<std::filesystem::path> kmc_path;

    std::filesystem::path kmc_father_path, kmc_mother_path{};


    format_fasta::file_extensions.push_back("fa");
    format_fasta::file_extensions.push_back("fasta");

    myparser.add_option(kmc_path, 'k', "kmc-db", "Unified KMC file of parents (can give multiple, shorter k-mers first).");

    // TODO: make sure at least one/two path(s) are passed.
    myparser.add_option(reads_child_path1, 'x', "pe1", "Please provide the read file of the child (1. PE)", option_spec::DEFAULT, input_file_validator{{"fa", "fasta", "fq","fastq","gz"}});
    myparser.add_option(reads_child_path2, 'y', "pe2", "Please provide the read file of the child (2. PE)", option_spec::DEFAULT, input_file_validator{{"fa", "fasta", "fq","fastq","gz"}});
    myparser.add_option(reads_child_path, 'z', "sp", "Please provide the read file of the child.", option_spec::DEFAULT, input_file_validator{{"fa", "fasta", "fq","fastq","gz"}});
    myparser.add_option(out_path, 'o', "out", "Please provide the output read file name.", option_spec::DEFAULT, output_file_validator{{"fq","fastq"}});

    // myparser.add_option(kmc_father_path, 'f', "kmc_father", "Please provide the KMC db of the father.", option_spec::DEFAULT);
    // myparser.add_option(kmc_mother_path, 'm', "kmc_mother", "Please provide the KMC db of the mother.", option_spec::DEFAULT);

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

    std::vector<CKMCFile> kmc_db;
    kmc_db.resize(kmc_path.size());
    for (uint32_t i = 0; i < kmc_path.size(); ++i)
    {
        if (!kmc_db[i].OpenForRA(kmc_path[i]))
        {
            std::cerr << "Could not open KMC file.\n";
            return 1;
        }
    }

    // CKMCFile kmc_father, kmc_mother;
    // if (!kmc_father.OpenForRA(kmc_father_path) || !kmc_mother.OpenForRA(kmc_mother_path))
    // {
    //     std::cerr << "Could not open KMC files of parents.\n";
    //     return 1;
    // }

    if (paired_end)
        run<paired_sequence_file_output>(reads_child_path1, reads_child_path2, kmc_db, /*kmc_father, kmc_mother,*/ out_path, threads);
    else
        run<single_sequence_file_output>(reads_child_path, reads_child_path, kmc_db, /*kmc_father, kmc_mother,*/ out_path, threads);

    std::cout << '\n';
    return 0;
}
