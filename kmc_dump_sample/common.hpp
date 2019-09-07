#pragma once

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

struct single_sequence_file_output {

    using fields_type = fields<field::SEQ, field::ID, field::QUAL>;
    using formats_type = type_list<format_fastq>; // format_fasta,
    using stream_char_type = char;

    sequence_file_output<fields_type, formats_type, stream_char_type> sf;

    single_sequence_file_output(std::filesystem::path && path): sf(path.string() + std::string{".fq"})
    {
        sf.options.fasta_letters_per_line = 0;
    }

    void push_back(auto && r, auto && /*read_unused*/)
    {
        sf.push_back(r);
    }
};

std::string get_perc(auto && count, auto && total)
{
    std::string perc{std::to_string(truncf((count / static_cast<float>(total))*10000)/100)};
    perc.erase(perc.find_last_not_of('0') + 1, std::string::npos);
    return perc;
}

enum class fo_mo_state : uint8_t { NONE = 0, FO_LAST = 1, MO_LAST = 2 };

template <typename t1, typename t2, typename t3, typename t4>
inline void loadReads(t1 & chunks_container1, t2 & chunks_container2, t3 & it1, t4 & it2)
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
                    chunks_container1.emplace_back(*chunk_it);
            }

            #pragma omp section
            {
                chunks_container2.clear();
                auto chunk = *it2;
                for (auto chunk_it = chunk.begin(); chunk_it != chunk.end(); chunk_it++)
                    chunks_container2.emplace_back(*chunk_it);
            }
        }
    }
}

template <typename t1, typename t2>
inline void loadReads(t1 & chunks_container1, t2 & it1)
{
    chunks_container1.clear();
    auto chunk = *it1;
    for (auto chunk_it = chunk.begin(); chunk_it != chunk.end(); chunk_it++)
        chunks_container1.emplace_back(*chunk_it);
}
