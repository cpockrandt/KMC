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
    std::filesystem::path read_path{};
    std::filesystem::path kmc_path{};

    int threads{omp_get_max_threads()};
};

inline void run(program_options const & options, CKMCFile & kmc_db)
{
    uint32 kmer_length;
    {
        uint32 _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
        uint64 _max_count, _total_kmers;
        kmc_db.Info(kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
    }

    uint64_t no_reads = 0;
    uint64_t no_kmers = 0;

    std::vector<std::vector<uint32>> counters(options.threads);
    for (int i = 0; i < options.threads; ++i)
    {
        counters[i].reserve(15'000);
    }

    sequence_file_input<sequence_file_input_default_traits_dna> fin{options.read_path};
    auto chunks = fin | ranges::views::chunk(100);
    std::vector<typename decltype(fin)::record_type> chunks_container;
    chunks_container.reserve(100);
    auto it{chunks.begin()};

    std::array<uint32_t, 255> occ{0};

    for (; it != chunks.end(); ++it)
    {
        loadReads(chunks_container, it);

        #pragma omp parallel for num_threads(options.threads)
        for (uint64_t id = 0; id < chunks_container.size(); ++id)
        {
            auto & read{chunks_container[id]};
            auto read_view = get<field::SEQ>(read) | views::to_char;
            std::string read_str(read_view.begin(), read_view.end());
            auto & counter = counters[omp_get_thread_num()];

            kmc_db.GetCountersForRead(read_str, counter);

            no_kmers += counter.size();

            for (uint32_t i = 0; i < counter.size(); ++i)
                occ[counter[i] - 1]++;
        }

        no_reads += chunks_container.size();

        if (no_reads > chunks_container.size())
           std::cout << "\033[" << (3 + 2*11) << "A";

        // for (uint64_t i = 0; i < 255; ++i)
        //     std::cout << (i + 1) << ' ';
        // std::cout << '\n';
        //
        // for (uint64_t i = 0; i < 255; ++i)
        //     std::cout << occ[i] << ' ';
        // std::cout << '\n';

        std::cout << "k-mers  1x:\t" << occ[1 - 1] << " (" << (100.0f * occ[1 - 1] / no_kmers) << " %)\n";
        std::cout << "k-mers  2x:\t" << occ[2 - 1] << " (" << (100.0f * occ[2 - 1] / no_kmers) << " %)\n";
        std::cout << "k-mers  3x:\t" << occ[3 - 1] << " (" << (100.0f * occ[3 - 1] / no_kmers) << " %)\n";

        for (uint32_t i = 10; i <= 20; ++i)
            std::cout << "k-mers " << i << "x:\t" << occ[i - 1] << " (" << (100.0f * occ[i - 1] / no_kmers) << " %)\n";
        for (uint32_t i = 25; i <= 35; ++i)
            std::cout << "k-mers " << i << "x:\t" << occ[i - 1] << " (" << (100.0f * occ[i - 1] / no_kmers) << " %)\n";

        std::cout.flush();
    }
}

int _tmain(int argc, char* argv[])
{
    argument_parser parser{"Heterozygosity-estimation", argc, argv};

    parser.info.author = "Christopher Pockrandt";
    parser.info.short_description = "Estimating the heterozygosity.";
    parser.info.version = "0.0.1";

    program_options options;
    parser.add_positional_option(options.read_path, "Read file.", input_file_validator{{"fa", "fasta", "fq","fastq","gz"}});
    parser.add_positional_option(options.kmc_path, "KMC file.");
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

    CKMCFile kmc_db;
    if (!kmc_db.OpenForRA(options.kmc_path))
    {
        std::cerr << "Could not open KMC file.\n";
        return 1;
    }

    std::cout << '\n' << std::fixed << std::setprecision(2);

    run(options, kmc_db);

    return 0;
}
