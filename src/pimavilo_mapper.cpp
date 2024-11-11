#include <iostream>
#include <getopt.h>
#include <cstdlib>

void print_help() {
    std::cout << "Usage: pimavilo_mapper [options] <file1> <file2>\n"
              << "Options:\n"
              << "  -h, --help       Show this help message and exit\n"
              << "  --version        Show version information\n";
}

void print_version() {
    std::cout << "pimavilo_mapper version " << PROJECT_VERSION << "\n";
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
                // getopt_long already prints an error message
                print_help();
                return 1;
            default:
                abort();
        }
    }

    // After options, we expect two additional arguments (file1 and file2)
    if (argc - optind < 2) {
        std::cerr << "Error: Missing required file arguments.\n";
        print_help();
        return 1;
    }

    // Read the file arguments
    const char* file1 = argv[optind];
    const char* file2 = argv[optind + 1];

    // Print file names (for demonstration purposes)
    std::cout << "Processing files: " << file1 << " and " << file2 << "\n";

    // Add your processing code here

    return 0;
}
