#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdlib>
#include <algorithm>
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include <cstddef>
#include <cstdlib>
#include "pimavilo/Align.h"
#include "pimavilo/Minimizers.hpp"
#include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <cmath>

#define EPS 700
#define CLUSTER 3

using namespace std;

// Default values for alignment scores
int match_score = 2;
int mismatch_penalty = -1;
int gap_penalty = -2;

// Set k, w and f parameters
unsigned int kmer_len = 5;
unsigned int window_len = 15;
double frequency_threshold = 0.001;

bool cigar = false;

// Default alignment type
pimavilo::AlignmentType alignment_type = pimavilo::AlignmentType::Local;

void print_help()
{
    cerr << "Usage: pimavilo_mapper [options] <reference_file> <fragments_file>\n"
         << "Options:\n"
         << "  -h, --help       Show this help message and exit\n"
         << "  --version        Show version information\n"
         << "  -m, --match      Value for match score (default: 2)\n"
         << "  -n, --mismatch   Value for mismatch penalty (default: -1)\n"
         << "  -g, --gap        Value for gap penalty (default: -2)\n"
         << "  -a, --alignment  Alignment type (can be: global, semiglobal, local; default: local)\n"
         << "  -k, --kmer_len   Kmer length (default: 5)\n"
         << "  -w, --window_len Window length (default: 15)\n"
         << "  -f, --frequency  Frequency threshold (default: 0.001)\n"
         << "  -c, --cigar      Use cigar string\n";
}

void print_version()
{
    cerr << "pimavilo_mapper version " << PROJECT_VERSION << "\n";
}

// Struct to represent a sequence and its length
struct Sequence_Fasta
{
private:
    const char *seq_name;
    const char *seq_data;
    std ::uint32_t name_len;
    std ::uint32_t data_len;

public:
    std::string name;
    std::string data;

    Sequence_Fasta(const char *seq_name, std::uint32_t name_len, const char *seq_data, std::uint32_t data_len) : name(seq_name, name_len), data(seq_data, data_len), data_len(data_len) {}
    std::size_t length() const { return data_len; }
};

struct Sequence_Fastq
{
private:
    const char *seq_name;
    const char *seq_data;
    const char *seq_qual;
    std ::uint32_t name_len;
    std ::uint32_t data_len;
    std ::uint32_t qual_len;

public:
    std::string name;
    std::string data;
    std::string qual;

    Sequence_Fastq(const char *seq_name, std::uint32_t name_len, const char *seq_data, std::uint32_t data_len, const char *seq_qual, std::uint32_t qual_len) : name(seq_name, name_len), data(seq_data, data_len), qual(seq_qual, qual_len), data_len(data_len) {}
    std::size_t length() const { return data_len; }
};

// Calculate statistics and print to stderr
template <typename Sequence>
void calculate_statistics(const std::vector<std::unique_ptr<Sequence>> &sequences)
{
    if (sequences.empty())
    {
        cerr << "No sequences found.\n";
        return;
    }

    std::size_t num_sequences = sequences.size();
    std::size_t total_length = 0;
    std::size_t min_length = sequences.begin()->get()->length();
    std::size_t max_length = sequences.begin()->get()->length();
    std::vector<size_t> lengths;

    for (const auto &seq : sequences)
    {
        std::size_t len = seq->length();
        total_length += len;
        min_length = min(min_length, len);
        max_length = max(max_length, len);
        lengths.push_back(len);
    }

    // Average length
    double average_length = static_cast<double>(total_length) / num_sequences;

    // N50 calculation
    std::sort(lengths.begin(), lengths.end(), std::greater<std::size_t>());
    std::size_t half_total_length = total_length / 2;
    std::size_t n50 = 0;
    std::size_t cumulative_length = 0;

    for (const auto &len : lengths)
    {
        cumulative_length += len;
        if (cumulative_length >= half_total_length)
        {
            n50 = len;
            break;
        }
    }

    // Output statistics to stderr
    cerr << "Number of sequences: " << num_sequences << "\n";
    cerr << "Total length of sequences: " << total_length << "\n";
    cerr << "Average sequence length: " << average_length << "\n";
    cerr << "Minimum sequence length: " << min_length << "\n";
    cerr << "Maximum sequence length: " << max_length << "\n";
    cerr << "N50 length: " << n50 << "\n";
}

int CountMsInCigar(const std::string &cigar) {
    int count = 0;
    size_t i = 0;

    while (i < cigar.size()) {
        int num = 0;

        // Parse the number preceding the operation
        while (i < cigar.size() && std::isdigit(cigar[i])) {
            num = num * 10 + (cigar[i] - '0');
            i++;
        }

        // Check if the operation is 'M'
        if (i < cigar.size() && cigar[i] == 'M') {
            count += num;  // Add the parsed number to the count
        }

        // Move to the next character
        i++;
    }

    return count;
}

//function to generate PAF output line
std::string GeneratePAF(const std::string &query_name, unsigned int query_len, string query,
                        unsigned int query_start, unsigned int query_end, char strand, 
                        const std::string &target_name, unsigned int target_len, string target,
                        unsigned int target_start, unsigned int target_end,
                        unsigned int alignment_length, unsigned int mapping_quality, const bool &cigar) {
        
    string cigar_string = "";
    unsigned int target_begin = 0;
    string str1, str2;
    if (query_start < 6 || query_end > query_len - 10) str1 = query;
    else str1 = query.substr(query_start - 5, query_end + 10 - query_start);
    if (target_start < 6 || target_end > target_len - 10) str2 = target;
    else str2 = target.substr(target_start - 5, target_end + 10 - target_start);
    int align = pimavilo::Align(str1.c_str(), str1.size(), 
                 str2.c_str(), str2.size(), alignment_type,
                 match_score, mismatch_penalty, gap_penalty, &cigar_string, &target_begin);

    int matching_bases = CountMsInCigar(cigar_string);

    std::string paf_line = query_name + "\t" + std::to_string(query_len) + "\t" +
                            std::to_string(query_start) + "\t" + std::to_string(query_end) + "\t" +
                            strand + "\t" + target_name + "\t" + std::to_string(target_len) + "\t" +
                            std::to_string(target_start) + "\t" + std::to_string(target_end) + "\t" +
                            std::to_string(matching_bases) + "\t" + std::to_string(alignment_length) + "\t" +
                            std::to_string(mapping_quality);

    if(cigar){
        paf_line += "\tcg:Z:" + cigar_string;
    }
    return paf_line;
}

bool areTheyClose(tuple<int,int,int> first_match, tuple<int,int,int> second_match, bool same_strand) { 
    if (same_strand) {
        return abs((get<1>(first_match) - get<2>(first_match))-(get<1>(second_match) - get<2>(second_match))) <= EPS;
    }
    else{
        return abs((get<1>(first_match) + get<2>(first_match))-(get<1>(second_match) + get<2>(second_match))) <= EPS;
    }
}

vector<int> LIS(vector<tuple<int,int,int>> matches){
    // Find LIS on reference indices while preserving fragment order
    vector<int> lis; // Stores indices of LIS elements
    vector<int> prev(matches.size(), -1); // For reconstructing the LIS

    for (size_t i = 0; i < matches.size(); ++i) {
        int ref_index = get<2>(matches[i]);

        // Find position in LIS using binary search
        auto pos = lower_bound(lis.begin(), lis.end(), i, [&](int lis_index, int val) {
            return get<2>(matches[lis_index]) < get<2>(matches[val]) && get<1>(matches[lis_index]) < get<1>(matches[val]);
        });

        size_t idx = distance(lis.begin(), pos);

        if (pos == lis.end()) {
            lis.push_back(i); // Extend the LIS
        } else {
            lis[idx] = i; // Replace to maintain smallest possible LIS
        }

        // Update previous element for LIS reconstruction
        if (idx > 0) {
            prev[i] = lis[idx - 1];
        }
    }

    // Reconstruct LIS from the indices
    vector<int> lis_indices;
    for (int i = lis.back(); i != -1; i = prev[i]) {
        lis_indices.push_back(i);
    }
    reverse(lis_indices.begin(), lis_indices.end());

    if(lis_indices.size() < 2){
        lis_indices.clear();
        return lis_indices;
    }
    // Output or process the LIS
    // if (!lis_indices.empty()) {
    //     const auto &first_match = matches[lis_indices.front()];
    //     const auto &last_match = matches[lis_indices.back()];
    // }

    return lis_indices;
}

struct TupleHash {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        auto h1 = std::hash<int>{}(std::get<0>(t));
        auto h2 = std::hash<int>{}(std::get<1>(t));
        auto h3 = std::hash<int>{}(std::get<2>(t));
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

int main(int argc, char *argv[]){
    // Define long options
    const struct option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"version", no_argument, nullptr, 'v'},
        {"match", required_argument, nullptr, 'm'},
        {"mismatch", required_argument, nullptr, 'n'},
        {"gap", required_argument, nullptr, 'g'},
        {"alignment_type", required_argument, nullptr, 'a'},
        {"kmer_len", required_argument, nullptr, 'k'},
        {"window_len", required_argument, nullptr, 'w'},
        {"frequency", required_argument, nullptr, 'f'},
        {"cigar", no_argument, nullptr, 'c'},
        {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    int opt;

    // Parse command line arguments
    while ((opt = getopt_long(argc, argv, "hvm:n:g:a:k:w:f:c", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'h':
            print_help();
            return 0;
        case 'v':
            print_version();
            return 0;
        case 'm':
            match_score = atoi(optarg);
            break;
        case 'n':
            mismatch_penalty = atoi(optarg);
            break;
        case 'g':
            gap_penalty = atoi(optarg);
            break;
        case 'a':
            if (strcmp(optarg, "global") == 0)
            {
                alignment_type = pimavilo::AlignmentType::Global;
            }
            else if (strcmp(optarg, "local") == 0)
            {
                alignment_type = pimavilo::AlignmentType::Local;
            }
            else if (strcmp(optarg, "semiglobal") == 0)
            {
                alignment_type = pimavilo::AlignmentType::SemiGlobal;
            }
            break;
        case 'k':
            kmer_len = atoi(optarg);
            break;
        case 'w':
            window_len = atoi(optarg);
            break;
        case 'f':
            frequency_threshold = atof(optarg);
            break;
        case 'c':
            cigar = true;
            break;
        case '?':
            print_help();
            return 1;
        default:
            abort();
        }
    }

    // After options, we expect two additional arguments (reference_file and fragments_file)
    if (argc - optind < 2)
    {
        cerr << "Error: Missing required file arguments.\n";
        print_help();
        return 1;
    }

    std::ios_base::sync_with_stdio(false);

    // Load sequences from reference and fragments files
    const std::string reference_file = argv[optind];
    const std::string fragments_file = argv[optind + 1];

    cerr << "Processing reference file: " << reference_file << "\n";
    cerr << "Processing fragments file: " << fragments_file << "\n";

    auto parser = bioparser::Parser<Sequence_Fasta>::Create<bioparser::FastaParser>(reference_file);
    auto reference_sequences = parser->Parse(-1); // parsing reference sequences
    cerr << "Name  of reference sequence: " << reference_sequences[0]->name << "\n";
    cerr << "Length of reference sequence: " << reference_sequences[0]->length() << "\n"; // Output of name and length of reference sequence
    cerr << "Parsing data...\n";
    parser = bioparser::Parser<Sequence_Fasta>::Create<bioparser::FastaParser>(fragments_file);
    auto sequences = parser->Parse(-1); // parsing fragment sequences

    auto frag_parser = bioparser::Parser<Sequence_Fasta>::Create<bioparser::FastaParser>(fragments_file);
    auto fragments_sequences = frag_parser->Parse(-1);

    cerr << "Creating reference minimizer index...\n";
    auto reference_minimizers = pimavilo::Minimize(reference_sequences[0]->data.c_str(), reference_sequences[0]->length(), kmer_len, window_len);
    auto reference_index = pimavilo::BuildMinimizerIndex(reference_minimizers, frequency_threshold);

    cerr << "Calculating frequencies...\n";
    // process the reference file
    std::unordered_map<unsigned int, unsigned int> minimizer_counts_ref;
    for (const auto &ref_seq : reference_sequences)
    {
        auto minimizers = pimavilo::Minimize(ref_seq->data.c_str(), ref_seq->length(), kmer_len, window_len);
        for (const auto &[hash, index, strand] : minimizers)
        {
            minimizer_counts_ref[hash]++;
        }
    }

    // process the fragments file
    std::unordered_map<unsigned int, unsigned int> minimizer_counts_frag;
    // omp_set_num_threads(16);
    // #pragma omp parallel for schedule(dynamic)
    cerr << "Mapping...\n";
    for (size_t frag_idx = 0; frag_idx < fragments_sequences.size(); ++frag_idx) {
        // cerr << "frag_idx " << frag_idx << "\n";
        const auto &frag_seq = fragments_sequences[frag_idx];
        auto minimizers = pimavilo::Minimize(frag_seq->data.c_str(), frag_seq->length(), kmer_len, window_len);

        // Find matches between fragment minimizers and reference index
        unordered_set<tuple<int, int, int>, TupleHash> lookup_same; // (hash, fragment_index, reference_index) in the same strand
        unordered_set<tuple<int, int, int>, TupleHash> lookup_diff; // (hash, fragment_index, reference_index) in different strands
        vector<tuple<int, int, int>> matches_same; // (hash, fragment_index, reference_index) in the same strand
        vector<tuple<int, int, int>> matches_diff; // (hash, fragment_index, reference_index) in different strands
        for (const auto &[hash, frag_index, strand] : minimizers) {
            minimizer_counts_frag[hash]++;
            if (reference_index.count(hash)) {
                for (const auto &[ref_index, ref_strand] : reference_index.at(hash)) {
                    if (ref_strand == strand && lookup_same.find(tuple<int, int, int>{hash, frag_index, ref_index}) == lookup_same.end()) {
                        lookup_same.emplace(tuple<int, int, int>{hash, frag_index, ref_index});
                        matches_same.emplace_back(hash, frag_index, ref_index);
                    }
                    else if(ref_strand != strand && lookup_diff.find(tuple<int, int, int>{hash, frag_index, ref_index}) == lookup_diff.end()) {
                        lookup_diff.emplace(tuple<int, int, int>{hash, frag_index, ref_index});
                        matches_diff.emplace_back(hash, frag_index, ref_index);
                    }
                }
            }
        }
        vector<tuple<int,int,int>> cluster;
        if (!matches_same.empty()) {
            cluster.emplace_back(matches_same[0]);
            for (int i = 1; i < matches_same.size(); ++i) {
                if (!cluster.empty() && areTheyClose(matches_same[i], cluster.back(), true)) {
                    cluster.emplace_back(matches_same[i]);
                }
                else{
                    if (cluster.size() >= CLUSTER){
                        vector<int> s = LIS(cluster);
                        if (!s.empty()){
                            cerr << GeneratePAF(frag_seq->name, frag_seq->length(), frag_seq->data, 
                            get<1>(cluster[s.front()]), get<1>(cluster[s.back()]), 'F',
                            reference_sequences[0]->name, reference_sequences[0]->length(), reference_sequences[0]->data,
                            get<2>(cluster[s.front()]), get<2>(cluster[s.back()]), 
                            max(get<1>(cluster[s.back()]) - get<1>(cluster[s.front()]), get<2>(cluster[s.back()]) - get<2>(cluster[s.front()])),
                            100, cigar) << "\n";
                        }
                    }
                    cluster.clear();
                    cluster.emplace_back(matches_same[i]);
                }
            }
        }
        if (cluster.size() >= CLUSTER) {
            vector<int> s = LIS(cluster);
            if (!s.empty()){
                cerr << GeneratePAF(frag_seq->name, frag_seq->length(), frag_seq->data, 
                get<1>(cluster[s.front()]), get<1>(cluster[s.back()]), 'F',
                reference_sequences[0]->name, reference_sequences[0]->length(), reference_sequences[0]->data,
                get<2>(cluster[s.front()]), get<2>(cluster[s.back()]), 
                max(get<1>(cluster[s.back()]) - get<1>(cluster[s.front()]), get<2>(cluster[s.back()]) - get<2>(cluster[s.front()])),
                100, cigar) << "\n";
            }
        }
        cluster.clear();

        if (matches_diff.empty()) continue;
        cluster.emplace_back(matches_diff[0]);
        for (int i = 1; i < matches_diff.size(); ++i) {
            if (!cluster.empty() && areTheyClose(matches_diff[i], cluster.back(), false)) {
                cluster.emplace_back(matches_diff[i]);
            }
            else{
                if (cluster.size() >= CLUSTER){
                    vector<int> s = LIS(cluster);
                    if (!s.empty()){
                        cerr << GeneratePAF(frag_seq->name, frag_seq->length(), frag_seq->data, 
                        get<1>(cluster[s.front()]), get<1>(cluster[s.back()]), 'F',
                        reference_sequences[0]->name, reference_sequences[0]->length(), reference_sequences[0]->data,
                        get<2>(cluster[s.front()]), get<2>(cluster[s.back()]), 
                        max(get<1>(cluster[s.back()]) - get<1>(cluster[s.front()]), get<2>(cluster[s.back()]) - get<2>(cluster[s.front()])),
                        100, cigar) << "\n";
                    }
                }
                cluster.clear();
                cluster.emplace_back(matches_diff[i]);
            }
        }
        if (cluster.size() >= CLUSTER) {
            vector<int> s = LIS(cluster);
            if (!s.empty()){
                cerr << GeneratePAF(frag_seq->name, frag_seq->length(), frag_seq->data, 
                get<1>(cluster[s.front()]), get<1>(cluster[s.back()]), 'F',
                reference_sequences[0]->name, reference_sequences[0]->length(), reference_sequences[0]->data,
                get<2>(cluster[s.front()]), get<2>(cluster[s.back()]), 
                max(get<1>(cluster[s.back()]) - get<1>(cluster[s.front()]), get<2>(cluster[s.back()]) - get<2>(cluster[s.front()])),
                100, cigar) << "\n";
            }
        }
    }
        
    // Anayze Minimizer Statistics
    auto analyze_minimizers = [](const std::unordered_map<unsigned int, unsigned int> &counts, double f)
    {
        size_t total_minimizers = counts.size();
        size_t singletons = std::count_if(counts.begin(), counts.end(), [](const std::pair<const unsigned int, unsigned int> &entry)
                                          { return entry.second == 1; });

        std::vector<unsigned int> frequencies;
        for (const auto &entry : counts)
        {
            frequencies.push_back(entry.second);
        }

        std::sort(frequencies.begin(), frequencies.end(), std::greater<unsigned int>());

        size_t exclude_count = static_cast<size_t>(f * frequencies.size());
        size_t max_freq = (frequencies.size() > exclude_count) ? frequencies[exclude_count] : 0;

        // print statistics
        cerr << "Number of minimizers: " << total_minimizers << "\n";
        cerr << "Fraction of singletons: " << singletons << "\n";
        cerr << "Max frequency: " << max_freq << "\n";
    };

    // Print statistics for both reference and fragment minimizers
    std::cerr << "Reference Minimizer Statistics:\n";
    analyze_minimizers(minimizer_counts_ref, frequency_threshold);

    std::cerr << "\nFragment Minimizer Statistics:\n";
    analyze_minimizers(minimizer_counts_frag, frequency_threshold);

    std::cerr << "\n";
    // Calculate and output statistics for fragment sequences
    calculate_statistics<Sequence_Fasta>(sequences);

    

    unsigned int target_begin = 0;
    string cigar_string = "";
    // Perform dummy alignment
    cerr << "\nDummy alignment " << pimavilo::Align(sequences[0]->data.c_str(), sequences[0]->length(), sequences[2]->data.c_str(), sequences[2]->length(), alignment_type, match_score, mismatch_penalty, gap_penalty, &cigar_string, &target_begin) << "\n";
    cerr << "begin: " << target_begin << " cigar: " << cigar_string << endl;

    //Generate PAF output
    cerr << "\nDummy PAF:\n";
    std::string paf_line = GeneratePAF(sequences[0]->name, sequences[0]->length(), sequences[0]->data, 
            1, sequences[0]->length(), 'F', sequences[2]->name, sequences[2]->length(), 
            sequences[2]->data, 1, sequences[2]->length(), 100, 100, &cigar_string);
    cerr << paf_line << endl;
    
    return 0;
}
