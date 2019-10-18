#include <iostream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/all.hpp>

#include <seqan3/range/views/to.hpp>

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

    uint64_t superreads_father_count = 0;
    uint64_t superreads_mother_count = 0;
    uint64_t superreads_bases_father_count = 0;
    uint64_t superreads_bases_mother_count = 0;

    std::vector<std::vector<dna4>> superreads;
    superreads.reserve(9'025'510 + 4'602'915);

    std::cout << mytime() << "Reading superreads of father ..." << std::endl;
    sequence_file_input fin_father{path_superreads_father/*, format_fasta{}*/};
    for (auto & rec : fin_father)
    {
        auto seq = get<field::SEQ>(rec)
                 | std::views::transform([](dna5 c){ return static_cast<dna4>(c); })
                 | views::to<std::vector<dna4> >;
        superreads.push_back(seq);
        superreads_bases_father_count += seq.size();
    }
    superreads_father_count += superreads.size();

    std::cout << mytime() << "Reading superreads of mother ..." << std::endl;
    sequence_file_input fin_mother{path_superreads_mother/*, format_fasta{}*/};
    for (auto & rec : fin_mother)
    {
        auto seq = get<field::SEQ>(rec)
                 | std::views::transform([](dna5 c){ return static_cast<dna4>(c); })
                 | views::to<std::vector<dna4> >;
        superreads.push_back(seq);
        superreads_bases_mother_count += seq.size();
    }
    superreads_mother_count = superreads.size() - superreads_father_count;

    std::cout << "- Bases: " << superreads_bases_father_count << ", " << superreads_bases_mother_count << std::endl;
    std::cout << "- Reads: " << superreads_father_count << ", " << superreads_mother_count << std::endl;

    uint64_t const bits_father = superreads_bases_father_count + superreads_father_count;
    uint64_t const bits_mother = superreads_bases_mother_count + superreads_mother_count;

    std::cout << mytime() << "Building FM index ..." << std::endl;
    bwt_binner bwt(superreads, bits_father, bits_mother);

    std::cout << mytime() << "Storing data on disk ..." << std::endl;
    {
        std::ofstream os{path_index_out / "index", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(bwt);
    }
    std::cout << mytime() << "Done!" << std::endl;

    // std::cout << mytime() << "Building FM index ..." << std::endl;
    // bi_fm_index<dna4, text_layout::collection> index{superreads};
    //
    // std::cout << "Storing FM index ..." << std::endl;
    // {
    //     std::ofstream os{path_index_out / "fm_index", std::ios::binary};
    //     cereal::BinaryOutputArchive oarchive{os};
    //     oarchive(index);
    // }
    //
    // std::cout << "Building bit vector ..." << std::endl;
    // // 0 is father, 1 is mother
    // // there will be superreads.size() many sentinel characters (for each superread)
    // // bit_vector_length.size() will cover the sentinels for each superread, that will be at the beginning of the BWT
    // // uint64_t const bit_vector_length = superreads_bases_father_count + superreads_bases_mother_count + superreads.size();
    // sdsl::bit_vector haplotype_indicator(superreads_bases_father_count + superreads_father_count);
    // haplotype_indicator.flip(); // set all bits to 1 (temporarily)
    // haplotype_indicator.resize(haplotype_indicator.size() + superreads_bases_mother_count + superreads_mother_count); // extend to accommodate bases for mother (initialized with 0)
    // haplotype_indicator.flip(); // father is now set to 0, mother to 1.
    // debug_stream << haplotype_indicator << std::endl;
    //
    // std::cout << "Building interleaved bit vector ..." << std::endl;
    // // sdsl::rank_support_v<> rb(&haplotype_indicator);
    // sdsl::bit_vector_il haplotype_indicator_il(haplotype_indicator);
    //
    // std::cout << "Storing interleaved bit vector ..." << std::endl;
    // {
    //     std::ofstream os{path_index_out / "indicator", std::ios::binary};
    //     cereal::BinaryOutputArchive oarchive{os};
    //     oarchive(haplotype_indicator_il);
    // }

    return 0;
}
