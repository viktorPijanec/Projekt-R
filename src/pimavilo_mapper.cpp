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

using namespace std;

void print_help()
{
    cout << "Usage: pimavilo_mapper [options] <reference_file> <fragments_file>\n"
         << "Options:\n"
         << "  -h, --help       Show this help message and exit\n"
         << "  --version        Show version information\n";
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
void calculate_statistics(const std::vector<Sequence> &sequences)
{
    if (sequences.empty())
    {
        cerr << "No sequences found.\n";
        return;
    }

    std::size_t num_sequences = sequences.size();
    std::size_t total_length = 0;
    std::size_t min_length = sequences[0].length();
    std::size_t max_length = sequences[0].length();
    std::vector<size_t> lengths;

    for (const auto &seq : sequences)
    {
        std::size_t len = seq.length();
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
        {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    int opt;

    // Parse command line arguments
    while ((opt = getopt_long(argc, argv, "h", long_options, &option_index)) != -1)
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
    auto reference_sequences = parser->Parse(-1);
    auto reference_sequence = std::move(reference_sequences[0]);
    cout
        << reference_sequence->name << "\n"
        << reference_sequence->data.substr(0, 200) << "\n";                                     // ispisano prvih 200 znakova
    parser = bioparser::Parser<Sequence_Fasta>::Create<bioparser::FastaParser>(fragments_file); // procitani fragmenti
    auto sequences = parser->Parse(-1);
    // for (const auto &seq : sequences)
    // {
    //     cout << seq->name << "\n"
    //          << seq->data.substr(0, 200) << "\n";
    // }
    // auto reference_sequences = <bioparser::FastqParser, Sequence_Fastq>(fragments_file);

    // Output names and lengths of sequences in the reference file
    cerr
        << "Reference sequences:\n";
    for (const auto &seq : sequences)
    {
        cerr << "- " << seq->name << " (" << seq->length() << " bp)\n";
    }

    // // Calculate and output statistics for fragment sequences
    // cerr << "\nFragment sequence statistics:\n";
    // calculate_statistics(fragment_sequences);

    return 0;
}
