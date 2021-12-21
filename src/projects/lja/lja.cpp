#include "pipeline.hpp"



std::string constructMessage() {
    std::stringstream ss;
    ss << "LJA: genome assembler for PacBio HiFi reads based on de Bruijn graph.\n";
    ss << "Usage: lja [options] -o <output-dir> --reads <reads_file> [--reads <reads_file2> ...]\n\n";
    ss << "Basic options:\n";
    ss << "  -o <file_name> (or --output-dir <file_name>)  Name of output folder. Resulting graph will be stored there.\n";
    ss << "  --reads <file_name>                           Name of file that contains reads in fasta or fastq format. This option can be used any number of times in the same command line. In this case reads from all specified files will be used as an input.\n";
    ss << "  -h (or --help)                                Print this help message.\n";
    ss << "\nAdvanced options:\n";
    ss << "  -t <int> (or --threads <int>)                 Number of threads. The default value is 16.\n";
    ss << "  -k <int>                                      Value of k used for initial error correction.\n";
    ss << "  -K <int>                                      Value of k used for final error correction and initialization of multiDBG.\n";
    ss << "  --diploid                                     Use this option for diploid genomes. By default LJA assumes that the genome is haploid or inbred.\n";
    return ss.str();
}

int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "threads=16", "k-mer-size=501", "window=2000", "K-mer-size=5001", "Window=500",
                     "cov-threshold=3", "rel-threshold=10", "Cov-threshold=3", "Rel-threshold=7",
                     "unique-threshold=40000", "dump", "dimer-compress=32,32,1", "restart-from=none", "load",
                     "alternative", "diploid", "debug", "help"},
                    {"reads", "paths", "ref"},
                    {"o=output-dir", "t=threads", "k=k-mer-size","w=window", "K=K-mer-size","W=Window", "h=help"},
                    constructMessage());
    parser.parseCL(argc, argv);
    if (parser.getCheck("help")) {
        std::cout << parser.message() << std::endl;
        return 0;
    }
    if (!parser.check().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parser.check() << "\n" << std::endl;
        std::cout << parser.message() << std::endl;
        return 1;
    }

    bool debug = parser.getCheck("debug");
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    logger.info() << "Hello! You are running La Jolla Assembler (LJA), a tool for genome assembly from PacBio HiFi reads\n";
    logging::logGit(logger, dir / "version.txt");
    bool diplod = parser.getCheck("diploid");
    std::string first_stage = parser.getValue("restart-from");
    bool skip = first_stage != "none";
    bool load = parser.getCheck("load");
    logger.info() << "LJA pipeline started" << std::endl;

    size_t threads = std::stoi(parser.getValue("threads"));

    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library paths = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paths"));
    io::Library ref_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("ref"));
    if(!io::CheckLibrary(lib + paths +ref_lib)) {
        exit(1);
    }
    multigraph::LJAPipeline pipeline (ref_lib);
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    size_t w = std::stoi(parser.getValue("window"));
    double threshold = std::stod(parser.getValue("cov-threshold"));
    double reliable_coverage = std::stod(parser.getValue("rel-threshold"));
    std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> corrected1;
    if (first_stage == "alternative")
        skip = false;
    corrected1 = pipeline.AlternativeCorrection(logger, dir / ("k" + itos(k)), lib, {}, paths, threads, k, w,
                                  threshold, reliable_coverage, false, false, skip, debug, load);
    if (first_stage == "alternative" || first_stage == "none")
        load = false;
    size_t K = std::stoi(parser.getValue("K-mer-size"));
    size_t W = std::stoi(parser.getValue("Window"));

    double Threshold = std::stod(parser.getValue("Cov-threshold"));
    double Reliable_coverage = std::stod(parser.getValue("Rel-threshold"));
    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold"));

    if(first_stage == "phase2")
        skip = false;
    std::experimental::filesystem::path executable(argv[0]);
    std::experimental::filesystem::path py_path = executable.parent_path() / "run_rr.py";
    logger.trace() << "py_path set to " << py_path.string() << std::endl;
    std::vector<std::experimental::filesystem::path> corrected2 =
            pipeline.SecondPhase(logger, dir / ("k" + itos(K)), {corrected1.first}, {corrected1.second}, paths,
                        threads, K, W, Threshold, Reliable_coverage, unique_threshold, py_path, diplod, skip, debug, load);
    if(first_stage == "phase2")
        load = false;
    if(first_stage == "polishing")
        skip = false;
    std::vector<std::experimental::filesystem::path> uncompressed_results =
            pipeline.PolishingPhase(logger, threads, dir/ "uncompressing", dir, corrected2[1],
                           corrected2[0],
                            lib, StringContig::max_dimer_size / 2, K, skip, debug);
    if(first_stage == "polishing")
        load = false;
    logger.info() << "Final homopolymer compressed and corrected reads can be found here: " << corrected2[0] << std::endl;
    logger.info() << "Final graph with homopolymer compressed edges can be found here: " << corrected2[1] << std::endl;
    logger.info() << "Final graph can be found here: " << uncompressed_results[1] << std::endl;
    logger.info() << "Final assembly can be found here: " << uncompressed_results[0] << std::endl;
    logger.info() << "LJA pipeline finished" << std::endl;
    return 0;
}
