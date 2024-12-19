#include "Minimizers.hpp"
#include <limits>
#include <algorithm>

namespace pimavilo{

   std::string GetReverseComplement(const std::string& sequence){
      std::string complete(sequence.size() - 1);
      for(size_t i = 0; i < sequence.size() - 1; i++){
         switch(sequence[sequence.size() - 1 - i]){
            case 'A':
               complement[i] = 'T';
               break;
            case 'T':
               complement[i] = 'A';
               break;
            case 'C':
               complement[i] = 'G';
               break;
            case 'G':
               complement[i] = 'C';
               break;
            default:
               complement[i] = 'N';
               break;
         }
      }

      return complement;
   }

   std::vector<unsigned int> HashKmers(const std:: string& sequence, unsigned int kmer_len){
      std::vector<unsigned int> hashes;
      if(sequence.size() < kmer_len){
         return hashes;
      }

      for(size_t i = 0; i < sequence.size() - 1; i++){
         std::string kmer = sequence.substr(i, kmer_len);
         unsigned int hash = std::hash<std::string>{}(kmer);
         hashes.push_back(hash);
      }

      return hashes;
   }

   std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len){
      std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;

      if(sequence_len < kmer_len || window_len < kmer_len){
         throw std::invalid_argument("Invalid kmer length or window length");
      }

      std::string seq(sequence, sequence_len);
      std::string rev_comp = GetReverseComplement(seq);

      auto fown_hashes = HashKmers(seq, kmer_len);
      auto rev_comp_hashes = HashKmers(rev_comp, kmer_len);

      for(size_t i = 0; i < fown_hashes.size() - window_len; i++){
         unsigned int min_hash = std::numeric_limits<unsigned int>::max();
         unsigned int min_index = 0;
         bool is_minimizer = true;

         for(size_t j = i; j < i +window_len; j++){
            if(fown_hashes[j].first < min_hash){
               min_hash = fown_hashes[j].first;
               min_index = fown_hashes[j].second;
               is_minimizer = true;
            }
            if(rev_comp_hashes[j].first < min_hash){
               min_hash = rev_comp_hashes[j].first;
               min_index = rev_comp_hashes[j].second;
               is_minimizer = false;
            }
         }

         minimizers.emplace_back(min_hash, min_index, is_minimizer);
      }

      return minimizers;
   }
}