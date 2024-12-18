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

      //dovrsit treba

      return minimizers;
   }
}