#include <iostream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/all.hpp>

#include <seqan3/std/ranges>

#include "bwt_binning_common.hpp"

using namespace seqan3;

int main(int argc, char* argv[])
{
    argument_parser parser{"Haplotype-binning-bwt-construction", argc, argv};

    parser.info.author = "Christopher Pockrandt";
    parser.info.short_description = "BWT construction for trio binning.";
    parser.info.version = "0.0.1";

    std::filesystem::path path_superreads_father{}, path_superreads_mother{}, path_index_out{};

    parser.add_positional_option(path_superreads_father, "Please provide the superread file of the father.", input_file_validator{{"fa", "fasta", "gz"}});
    parser.add_positional_option(path_superreads_mother, "Please provide the superread file of the mother.", input_file_validator{{"fa", "fasta", "gz"}});
    parser.add_positional_option(path_index_out, "Path where the index is written to.");

    try
    {
        parser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    std::vector<std::vector<dna4>> superreads_f, superreads_m;
    superreads_f.reserve(9'025'510);
    superreads_m.reserve(4'602'915);

    std::cout << mytime() << "Reading superreads of father ..." << std::endl;
    sequence_file_input fin_father{path_superreads_father/*, format_fasta{}*/};
    for (auto & rec : fin_father)
    {
        auto seq = get<field::SEQ>(rec) | std::views::transform([](dna5 c){ return static_cast<dna4>(c); });
        superreads_f.push_back(seq);
    }

    std::cout << mytime() << "Reading superreads of mother ..." << std::endl;
    sequence_file_input fin_mother{path_superreads_mother/*, format_fasta{}*/};
    for (auto & rec : fin_mother)
    {
        auto seq = get<field::SEQ>(rec) | std::views::transform([](dna5 c){ return static_cast<dna4>(c); });
        superreads_m.push_back(seq);
    }

    std::cout << mytime() << "Building FM index of father ..." << std::endl;
    bi_fm_index<dna4, text_layout::collection, my_sdsl_wt_index_type> index_f(superreads_f);

    std::cout << mytime() << "Storing data on disk ..." << std::endl;
    {
        std::ofstream os{path_index_out / "index_f", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index_f);
    }


    std::cout << mytime() << "Building FM index of mother ..." << std::endl;
    bi_fm_index<dna4, text_layout::collection, my_sdsl_wt_index_type> index_m(superreads_m);

    std::cout << mytime() << "Storing data on disk ..." << std::endl;
    {
        std::ofstream os{path_index_out / "index_m", std::ios::binary};
	cereal::BinaryOutputArchive oarchive{os};
	oarchive(index_m);
    }
    std::cout << mytime() << "Done!" << std::endl;

    return 0;
}
