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

using namespace std;

void print_help()
{
    cout << "Usage: pimavilo_mapper [options] <reference_file> <fragments_file>\n"
         << "Options:\n"
         << "  -h, --help       Show this help message and exit\n"
         << "  --version        Show version information\n"
         << "  -m, --match      Value for match score (default: 2)\n"
         << "  -n, --mismatch   Value for mismatch penalty (default: -1)\n"
         << "  -g, --gap        Value for gap penalty (default: -2)\n"
         << "  -a, --alignment  Alignment type (can be: global, semiglobal, local; default: local)\n";
}

void print_version()
{
    cout << "pimavilo_mapper version " << PROJECT_VERSION << "\n";
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

int main(int argc, char *argv[])
{
    // Define long options
    const struct option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"version", no_argument, nullptr, 'v'},
        {"match", required_argument, nullptr, 'm'},
        {"mismatch", required_argument, nullptr, 'n'},
        {"gap", required_argument, nullptr, 'g'},
        {"alignment_type", required_argument, nullptr, 'a'},
        {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    int opt;

    // Default values for alignment scores
    int match_score = 2;
    int mismatch_penalty = -1;
    int gap_penalty = -2;

    // Default alignment type
    pimavilo::AlignmentType alignment_type = pimavilo::AlignmentType::Local;

    // Parse command line arguments
    while ((opt = getopt_long(argc, argv, "hvm:n:g:a:", long_options, &option_index)) != -1)
    {
        // cout << opt << endl;
        switch (opt)
        {
        case 'h':
            print_help();
            return 0;
        case 'v':
            print_version();
            return 0;
        case 'm':
            match_score = std::atoi(optarg);
            break;
        case 'n':
            mismatch_penalty = std::atoi(optarg);
            break;
        case 'g':
            gap_penalty = std::atoi(optarg);
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

    // Load sequences from reference and fragments files
    const std::string reference_file = argv[optind];
    const std::string fragments_file = argv[optind + 1];

    cerr << "Processing reference file: " << reference_file << "\n";
    cerr << "Processing fragments file: " << fragments_file << "\n";

    auto parser = bioparser::Parser<Sequence_Fasta>::Create<bioparser::FastaParser>(reference_file);
    auto reference_sequences = parser->Parse(-1); // parsing reference sequences
    auto reference_sequence = std::move(reference_sequences[0]);
    cerr << "Name  of reference sequence: " << reference_sequence->name << "\n";
    cerr << "Length of reference sequence: " << reference_sequence->length() << "\n"; // Output of name and length of reference sequence
    parser = bioparser::Parser<Sequence_Fasta>::Create<bioparser::FastaParser>(fragments_file);
    auto sequences = parser->Parse(-1); // parsing fragment sequences

    // Calculate and output statistics for fragment sequences
    calculate_statistics<Sequence_Fasta>(sequences);

    unsigned int target_begin = 0;
    string cigar_string = "";
    // Perform dummy alignment
    cerr << "Align " << pimavilo::Align(sequences[0]->data.c_str(), sequences[0]->length(), sequences[2]->data.c_str(), sequences[2]->length(), alignment_type, match_score, mismatch_penalty, gap_penalty, &cigar_string, &target_begin) << "\n";
    cerr << "begin: " << target_begin << " cigar: " << cigar_string << endl;
    return 0;
}
