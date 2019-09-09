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

inline void comp_win(auto const & qual, uint16_t const from, uint16_t const to, uint16_t & window_m, uint16_t & window_f)
{
    window_m = 0;
    window_f = 0;
    for (uint16_t i = from; i < to + 1; ++i)
    {
        if (qual[i] == '2'_phred42)
            ++window_m;
        else if (qual[i] == '1'_phred42)
            ++window_f;
    }
}

template <bool output_details>
inline void analyze_read_sliding_win(auto & read, uint16_t & mo, uint16_t & fo, uint16_t & fo_mo_switches)
{
    auto & qual = get<field::QUAL>(read);

    mo = 0;
    fo = 0;
    fo_mo_switches = 0;
    fo_mo_state last_fo_mo{fo_mo_state::NONE};

    uint16_t window_m{0};
    uint16_t window_f{0};
    comp_win(qual, 0, 30, window_m, window_f);
    // debug_stream << "init: " << window_f << " ... " << window_m << " ... " << std::endl;

    for (uint32_t i = 0; i < qual.size() - 31 + 2; ++i)
    {
        if (window_m >= 0.8f * 31 || window_f >= 0.8f * 31)
        {
            if (window_m > window_f)
            {
                ++mo;
                if (last_fo_mo == fo_mo_state::FO_LAST)
                    ++fo_mo_switches;
                last_fo_mo = fo_mo_state::MO_LAST;
            }
            else
            {
                ++fo;
                if (last_fo_mo == fo_mo_state::MO_LAST)
                    ++fo_mo_switches;
                last_fo_mo = fo_mo_state::FO_LAST;
            }
            i += 31 - 1; // i will be incremented by for loop!
            comp_win(qual, i + 1, i + 1 + 31 - 1, window_m, window_f);
            // debug_stream << "jump: " << window_f << " ... " << window_m << " ... " << std::endl;
            continue;
        }

        if (i < qual.size() - 31 + 1)
        {
            // if (i == 864)
            //     std::cout << "XXXXXXXXXXXXXx" << std::endl;

            // sliding window: increment if it is moving into the window
            if (qual[i + 31] == '2'_phred42)
                ++window_m;
            else if (qual[i + 31] == '1'_phred42)
                ++window_f;

            // sliding window: decrement if it is moving out of the window
            if (qual[i] == '2'_phred42)
                --window_m;
            else if (qual[i] == '1'_phred42)
                --window_f;
            // debug_stream << "slide (" << i << "): " << window_f << " ... " << window_m << " ... " << "( " << qual[i + 31] << " " << qual[i] << " )" << std::endl;
        }
    }

    // debug_stream << qual << std::endl << fo << " ... " << mo << " ... " << fo_mo_switches << std::endl;
    // exit(0);
}
