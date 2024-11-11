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

void print_help() {
    cout << "Usage: pimavilo_mapper [options] <reference_file> <fragments_file>\n"
         << "Options:\n"
         << "  -h, --help       Show this help message and exit\n"
         << "  --version        Show version information\n";
}

void print_version() {
    cout << "pimavilo_mapper version " << PROJECT_VERSION << "\n";
}

// Struct to represent a sequence and its length
struct Sequence {
    std::string name;
    std::string data;

    Sequence(const std::string& name, const std::string& data) : name(name), data(data) {}
    std::size_t length() const { return data.size(); }
};

// Load sequences from a FASTA or FASTQ file and return them as a vector of Sequence objects
template <typename Parser>
std::vector<Sequence> load_sequences(const std::string& file_path) {
    std::vector<Sequence> sequences;
    auto parser = Parser::Create<Sequence>(file_path);
    parser->Parse([&](std::unique_ptr<Sequence> s) {
        sequences.emplace_back(move(*s));
    });
    return sequences;
}

// Calculate statistics and print to stderr
void calculate_statistics(const std::vector<Sequence>& sequences) {
    if (sequences.empty()) {
        cerr << "No sequences found.\n";
        return;
    }

    std::size_t num_sequences = sequences.size();
    std::size_t total_length = 0;
    std::size_t min_length = sequences[0].length();
    std::size_t max_length = sequences[0].length();
    std::vector<size_t> lengths;

    for (const auto& seq : sequences) {
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

    for (const auto& len : lengths) {
        cumulative_length += len;
        if (cumulative_length >= half_total_length) {
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

int main(int argc, char *argv[]) {
    // Define long options
    const struct option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"version", no_argument, nullptr, 'v'},
        {nullptr, 0, nullptr, 0}
    };

    int option_index = 0;
    int opt;

    // Parse command line arguments
    while ((opt = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
        switch (opt) {
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
    if (argc - optind < 2) {
        cerr << "Error: Missing required file arguments.\n";
        print_help();
        return 1;
    }

    // Load sequences from reference and fragments files
    const std::string reference_file = argv[optind];
    const std::string fragments_file = argv[optind + 1];

    cerr << "Processing reference file: " << reference_file << "\n";
    cerr << "Processing fragments file: " << fragments_file << "\n";

    auto reference_sequences = load_sequences<bioparser::FastaParser<Sequence>>(reference_file);
    auto fragment_sequences = load_sequences<bioparser::FastqParser<Sequence>>(fragments_file);

    // Output names and lengths of sequences in the reference file
    cerr << "Reference sequences:\n";
    for (const auto& seq : reference_sequences) {
        cerr << "- " << seq.name << " (" << seq.length() << " bp)\n";
    }

    // Calculate and output statistics for fragment sequences
    cerr << "\nFragment sequence statistics:\n";
    calculate_statistics(fragment_sequences);

    return 0;
}
