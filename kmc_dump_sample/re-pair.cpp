#include "omp.h"
#include "stdafx.h"
#include <iostream>
#include <iomanip>

#include <range/v3/view/take_while.hpp>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

int main(int argc, char* argv[])
{
    argument_parser myparser{"Re-pairing of paired-end reads", argc, argv};

    myparser.info.author = "Christopher Pockrandt";
    myparser.info.short_description = "Given two fastq files it outputs two paired files for reads for which both ends exists. Input reads can be in any order.";
    myparser.info.version = "0.0.1";

    std::filesystem::path pe1_path{}, pe2_path{};

    myparser.add_positional_option(pe1_path, "First unordered, incomplete PE file.");
    myparser.add_positional_option(pe2_path, "Second unordered, incomplete PE file.");

    try
    {
        myparser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    using record_t = typename sequence_file_input<>::record_type;
    sequence_file_input fin1{pe1_path}, fin2{pe2_path};
    std::vector<record_t> pe1, pe2;

	struct {
		bool operator()(record_t & a, record_t & b) const
		{
			return get<field::ID>(a) < get<field::ID>(b);
		}
	} custom_less;

    // #pragma omp parallel num_threads(2)
    {
        // #pragma omp sections
        {
            // #pragma omp section
            {
                for (auto & rec : fin1)
                    pe1.emplace_back(rec);
                std::cout << "Finished reading PE1." << std::endl;
                std::sort(pe1.begin(), pe1.end(), custom_less);
                std::cout << "Finished sorting PE1." << std::endl;
            }

            // #pragma omp section
            {
                for (auto & rec : fin2)
                    pe2.emplace_back(rec);
                std::cout << "Finished reading PE2." << std::endl;
                std::sort(pe2.begin(), pe2.end(), custom_less);
                std::cout << "Finished sorting PE2." << std::endl;
            }
        }
    }

    std::cout << "PE1: " << pe1.size() << '\n';
    std::cout << "PE2: " << pe2.size() << '\n';

    // debug_stream << pe1[0] << '\n';
    // debug_stream << pe2[0] << '\n' << '\n';
    //
    // debug_stream << pe1[1] << '\n';
    // debug_stream << pe2[1] << '\n' << '\n';
    //
    // debug_stream << pe1[2] << '\n';
    // debug_stream << pe2[2] << '\n' << '\n';

    std::string pe1_path_out{pe1_path}, pe2_path_out{pe2_path};

    pe1_path_out.insert(static_cast<std::string>(pe1_path_out).find_first_of('.'), ".PAIRED");
    pe2_path_out.insert(static_cast<std::string>(pe2_path_out).find_first_of('.'), ".PAIRED");

    sequence_file_output fout1{pe1_path_out}, fout2{pe2_path_out};

    uint64_t i = 0, j = 0, total = 0;
	while (true)
	{
        // auto key1 = get<field::ID>(pe1[i]) | ranges::view::take_while([] (char c) { return c != ' '; });
        // auto key2 = get<field::ID>(pe2[j]) | ranges::view::take_while([] (char c) { return c != ' '; });
        // {pe1.substr pe1[i].find_first_of(' ')

        auto & id1{get<field::ID>(pe1[i])};
        auto & id2{get<field::ID>(pe2[j])};

        std::string key1{id1.substr(0, id1.find_first_of(' '))};
        std::string key2{id2.substr(0, id2.find_first_of(' '))};

        // std::cout << key1 << " ... " << key2 << '\n';

		if (key1 == key2)
		// if (ranges::equal(key1, key2))
		{
            fout1.push_back(pe1[i]);
            fout2.push_back(pe2[j]);

			++i, ++j, ++total;
		}
        // else if (ranges::less(key1, key2))
		else if (key1 > key2)
			++j;
		else
			++i;

		if (i == pe1.size() || j == pe2.size())
			break;
	}

    std::cout << "Joint: " << total << '\n';

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
