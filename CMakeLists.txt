cmake_minimum_required (VERSION 3.2)
project (kmc_dump_sample CXX)

# set(CMAKE_BUILD_TYPE Release)

# Find SeqAn3
find_package (SeqAn3 REQUIRED
              HINTS "seqan3/build_system")
find_package (OpenMP REQUIRED)

# set (CXX11_FOUND "1")
# find_package (ZLIB)
# set (CMAKE_PREFIX_PATH "seqan/util/cmake")
# set (SEQAN_INCLUDE_PATH "seqan/include")
#find_package (OpenMP REQUIRED)
# find_package (SeqAn REQUIRED)
# include_directories (${SEQAN_INCLUDE_DIRS})
# add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your project.
add_executable (kmc_dump_sample kmc_dump_sample/kmc_dump_sample.cpp kmc_dump/nc_utils.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp)
add_executable (kmc_unique kmc_dump_sample/kmc_find_unique.cpp kmc_dump/nc_utils.cpp kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -m64 -static-libgcc -static-libstdc++ -fopenmp -pthread -std=c++14 ${SEQAN_CXX_FLAGS}") #  -fsanitize=undefined

# target_link_libraries (kmc_dump_sample ${SEQAN_LIBRARIES})
target_link_libraries (kmc_dump_sample seqan3::seqan3)
target_link_libraries (kmc_dump_sample pthread)
target_link_libraries (kmc_unique seqan3::seqan3)
target_link_libraries (kmc_unique pthread)
