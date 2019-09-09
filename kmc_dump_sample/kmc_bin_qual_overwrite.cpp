#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include "../kmc_api/kmc_file.h"

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>

#include "common.hpp"

using namespace seqan3;

inline void write_qualities(auto & read, auto & counter, auto & kmc_db, uint32 const kmer_length)
{
    auto read_view = get<field::SEQ>(read) | view::to_char;
    auto & qual = get<field::QUAL>(read);
    std::string read_str(read_view.begin(), read_view.end());
    kmc_db.GetCountersForRead(read_str, counter);

    qual.resize(read_str.size());
    for (uint32_t i = 0; i < read_str.size() - kmer_length + 1; ++i)
    {
        if (counter[i] == 2)
            qual[i].assign_char('2');
        else if (counter[i] == 1)
            qual[i].assign_char('1');
        else
            qual[i].assign_char('.');
    }
    for (uint32_t i = read_str.size() - kmer_length + 1; i < read_str.size(); ++i)
        qual[i].assign_char('.');
}

template <typename output_t>
inline void run(std::filesystem::path const & p1, std::filesystem::path & p2,
                auto & kmc_db, auto & kmc_father, auto & kmc_mother,
                int const threads)
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
    std::string p_out{p1};
    p_out.insert(static_cast<std::string>(p1).find_last_of('.'), ".qual" + std::to_string(kmer_length));
    output_t fout_all{p_out};

    uint64_t no_reads = 0;

    std::vector<std::vector<uint32>> counters(threads);
    for (int i = 0; i < threads; ++i)
    {
        counters[i].reserve(250 - kmer_length + 1);
    }

    sequence_file_input fin1{p1};
    auto chunks1 = fin1 | ranges::view::chunk(10'000);
    std::vector<typename sequence_file_input<>::record_type> chunks_container1;
    chunks_container1.reserve(10'000);
    auto it1{chunks1.begin()};

    if constexpr (!paired_end)
        p2 = p1;
    sequence_file_input fin2{p2};
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

        // #pragma omp parallel for num_threads(threads)
        for (uint64_t id = 0; id < chunks_container1.size(); ++id)
        {
            auto & read1{chunks_container1[id]};
            auto & read2{chunks_container2[id]};
            auto & counter = counters[omp_get_thread_num()];

            write_qualities(read1, counter, kmc_db, kmer_length);
            if constexpr (paired_end)
                write_qualities(read2, counter, kmc_db, kmer_length);

            #pragma omp critical(common)
            {
                fout_all.push_back(read1, read2);
            }
        }

        no_reads += chunks_container1.size();
        if (no_reads > chunks_container1.size())
            std::cout << "\033[1A";

        std::cout << "\nTotal            :\t" << no_reads << '\n';

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

    int threads{omp_get_max_threads()};
    std::filesystem::path kmc_path{}, reads_child_path{}, reads_child_path1{}, reads_child_path2{};

    std::filesystem::path kmc_father_path, kmc_mother_path{};

    myparser.add_positional_option(kmc_path, "Please provide the unified KMC file of the parents.");

    // TODO: make sure at least one/two path(s) are passed.
    myparser.add_option(reads_child_path1, 'x', "pe1", "Please provide the read file of the child (1. PE)", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});
    myparser.add_option(reads_child_path2, 'y', "pe2", "Please provide the read file of the child (2. PE)", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});
    myparser.add_option(reads_child_path, 'z', "sp", "Please provide the read file of the child.", option_spec::DEFAULT, input_file_validator{{"fq","fastq","gz"}});

    myparser.add_option(kmc_father_path, 'f', "kmc_father", "Please provide the KMC db of the father.", option_spec::DEFAULT);
    myparser.add_option(kmc_mother_path, 'm', "kmc_mother", "Please provide the KMC db of the mother.", option_spec::DEFAULT);

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
        run<paired_sequence_file_output>(reads_child_path1, reads_child_path2, kmc_db, kmc_father, kmc_mother, threads);
    else
        run<single_sequence_file_output>(reads_child_path, reads_child_path, kmc_db, kmc_father, kmc_mother, threads);

    std::cout << '\n';
    return 0;
}
