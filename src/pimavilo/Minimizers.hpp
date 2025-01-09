#ifndef PIMAVILO_MINIMIZERS_HPP
#define PIMAVILO_MINIMIZERS_HPP

#include <vector>
#include <tuple>
#include <string>

namespace pimavilo{
   //compute reverse complement for DNA/RNA
   std::string GetReverseComplement(const std::string& sequence);

   std::vector<unsigned int> HashKmers(const std::string& sequence, unsigned int kmer_len);

   //find minimizers
   std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len);

    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> BuildMinimizerIndex(
        const std::vector<std::tuple<unsigned int, unsigned int, bool>>& minimizers
        double frequency_treshold
    );
}

#endif