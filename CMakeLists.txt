cmake_minimum_required (VERSION 3.2)
project (kmc_dump_sample CXX)

# set(CMAKE_BUILD_TYPE Release)

# Find SeqAn3
find_package (SeqAn3 REQUIRED
              HINTS "seqan3/build_system")
find_package (OpenMP REQUIRED)

#set (CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

# set (CXX11_FOUND "1")
# find_package (ZLIB)
# set (CMAKE_PREFIX_PATH "seqan/util/cmake")
# set (SEQAN_INCLUDE_PATH "seqan/include")
#find_package (OpenMP REQUIRED)
# find_package (SeqAn REQUIRED)
# include_directories (${SEQAN_INCLUDE_DIRS})
# add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your project.
# link_directories(kmc_tools/libs)

add_executable (kmc_tools kmc_api/mmer.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp
                          kmc_tools/kmc_header.cpp kmc_tools/kmc_tools.cpp kmc_tools/nc_utils.cpp kmc_tools/parameters_parser.cpp kmc_tools/parser.cpp kmc_tools/tokenizer.cpp kmc_tools/fastq_filter.cpp kmc_tools/fastq_reader.cpp kmc_tools/fastq_writer.cpp kmc_tools/percent_progress.cpp)
add_executable (kmc_dump_sample kmc_dump_sample/kmc_dump_sample.cpp kmc_dump/nc_utils.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp)
add_executable (kmc_dump kmc_dump/kmc_dump.cpp kmc_dump/nc_utils.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp)
add_executable (kmc_unique kmc_dump_sample/kmc_find_unique.cpp kmc_dump/nc_utils.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp)
add_executable (kmc_cmp_reads kmc_dump_sample/kmc_cmp_reads.cpp kmc_dump/nc_utils.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp)
add_executable (kmc_keep_unique kmc_dump_sample/kmc_keep_unique.cpp)
add_executable (kmc_remove_Ns kmc_dump_sample/kmc_remove_Ns.cpp)
#add_executable (kmc kmer_counter/kmer_counter.cpp kmer_counter/mmer.cpp kmer_counter/mem_disk_file.cpp kmer_counter/rev_byte.cpp kmer_counter/bkb_writer.cpp kmer_counter/cpu_info.cpp kmer_counter/bkb_reader.cpp kmer_counter/fastq_reader.cpp kmer_counter/timer.cpp kmer_counter/develop.cpp kmer_counter/kb_completer.cpp kmer_counter/kb_storer.cpp kmer_counter/kmer.cpp kmer_counter/splitter.cpp kmer_counter/kb_collector.cpp kmer_counter/raduls_sse2.cpp kmer_counter/raduls_sse41.cpp kmer_counter/raduls_avx2.cpp kmer_counter/raduls_avx.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14 -lz -lbz2 -fopenmp ${SEQAN_CXX_FLAGS}") #  -fsanitize=undefined

# target_link_libraries (kmc_dump_sample ${SEQAN_LIBRARIES})
#target_link_libraries (kmc_tools seqan3::seqan3)
#target_link_libraries (kmc_tools pthread)
#target_link_libraries (kmc z)
#target_link_libraries (kmc bz2)
target_link_libraries (kmc_tools z)
target_link_libraries (kmc_tools bz2)
target_link_libraries (kmc_dump_sample seqan3::seqan3)
target_link_libraries (kmc_dump_sample pthread)
target_link_libraries (kmc_unique seqan3::seqan3)
target_link_libraries (kmc_unique pthread)
target_link_libraries (kmc_cmp_reads seqan3::seqan3)
target_link_libraries (kmc_cmp_reads pthread)
target_link_libraries (kmc_keep_unique seqan3::seqan3)
target_link_libraries (kmc_keep_unique pthread)
target_link_libraries (kmc_remove_Ns seqan3::seqan3)
target_link_libraries (kmc_remove_Ns pthread)
