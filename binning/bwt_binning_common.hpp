#pragma once

#include <seqan3/search/fm_index/all.hpp>
#include <sdsl/bit_vector_il.hpp>

using namespace seqan3;

std::string mytime()
{
    auto r = time(nullptr);
    auto c = ctime(&r);
    std::string buf(c);
    buf.insert(0, "[");
    buf.append("] ");
    buf.erase(remove(buf.begin(), buf.end(), '\n'), buf.end());
    return buf;
}

using my_sdsl_wt_index_type =
    sdsl::csa_wt<sdsl::wt_blcd<sdsl::bit_vector,
                               sdsl::rank_support_v<>,
                               sdsl::select_support_scan<>,
                               sdsl::select_support_scan<0>>,
                 4,
                 100000000,
                 sdsl::sa_order_sa_sampling<>,
                 sdsl::isa_sampling<>,
                 sdsl::plain_byte_alphabet>;

struct bwt_binner
{
    bi_fm_index<dna4, text_layout::collection, my_sdsl_wt_index_type> index;
    sdsl::bit_vector_il<> indicator;

    bwt_binner() {}

    bwt_binner(std::vector<std::vector<dna4>> const & superreads,
               uint64_t const bits_father,
               uint64_t const bits_mother) : index(superreads)
    {
        // sdsl::bit_vector haplotype_indicator(bits_father); // sets bits for father to 0
        // haplotype_indicator.flip(); // set all bits to 1 (temporarily)
        // haplotype_indicator.resize(bits_father + bits_mother); // extend to accommodate also bases for mother (initialized with 0)
        // haplotype_indicator.flip(); // father is now set to 0, mother to 1.
        // debug_stream << haplotype_indicator << std::endl;

        // auto cursor{index.begin()};
        // debug_stream << "range: " << cursor.fwd_lb << ", " << cursor.fwd_rb << std::endl;
        // debug_stream << "count: " << cursor.count() << std::endl;

        // std::cout << "- bits father: " << bits_father << std::endl;
        // std::cout << "- bits mother: " << bits_mother << std::endl;

        // fwd: mother comes first (since text is reversed)! [on rev it is reversed twice!]
        // for (uint64_t i = 0; i < index.fwd_fm.index.size(); ++i)
        //     std::cout << "fwd: " << index.fwd_fm.index[i] << std::endl;
        std::cout << mytime() << "Building bit vector ..." << std::endl;
        sdsl::bit_vector haplotype_indicator(bits_father + bits_mother);
        std::cout << std::fixed << std::setprecision(2);
        for (uint64_t i = 0; i < index.fwd_fm.index.size(); ++i)
        {
            if (index.fwd_fm.index[i] >= bits_mother)
                haplotype_indicator[i] = 0;
            else
                haplotype_indicator[i] = 1;
	    if (i % 4096 == 0)
	    {
	    	std::cout << "\rProgress: " << (100.0f * i / index.fwd_fm.index.size());
		std::cout.flush();
	    }
        }
        // debug_stream << haplotype_indicator << std::endl << std::endl;

        // for (uint64_t i = 0; i < index.rev_fm.index.size(); ++i)
        //     std::cout << "rev: " << index.rev_fm.index[i] << std::endl;
        // for (uint64_t i = 0; i < index.rev_fm.index.size(); ++i)
        // {
        //     if (index.rev_fm.index[i] >= bits_father)
        //         haplotype_indicator[i] = 1;
        //     else
        //         haplotype_indicator[i] = 0;
        // }
        // debug_stream << haplotype_indicator << std::endl << std::endl;

        // TODO: Set sampling rate to 1

        std::cout << mytime() << "Building interleaved bit vector ..." << std::endl;
        sdsl::bit_vector_il indicator_tmp(haplotype_indicator);
        std::swap(indicator, indicator_tmp);
    }

    template <class archive_t>
    void serialize(archive_t & archive)
    {
        archive(index);
        archive(indicator);
        //for (uint32_t i = 0; i < indicator.size(); ++i)
        //    std::cout << indicator[i] << std::endl;
    }
};
