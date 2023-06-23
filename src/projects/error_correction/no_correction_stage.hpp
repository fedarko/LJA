#pragma once

std::unordered_map<std::string, std::experimental::filesystem::path>
NoCorrection(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const io::Library &reads_lib,
             const io::Library &pseudo_reads_lib, const io::Library &paths_lib, size_t k, size_t w, bool debug, bool load) {
    logger.info() << "Performing graph construction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    hashing::RollingHash hasher(k);
    ensure_dir_existance(dir);
    io::Library construction_lib = reads_lib + pseudo_reads_lib;
    SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, construction_lib, dir, threads, (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string()) :
                    DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
    dbg.fillAnchors(w, logger, threads);
    size_t extension_size = std::max<size_t>(k * 2, 1000);
    ReadLogger readLogger(threads, dir/"read_log.txt");
    RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, true, false);
    RecordStorage extra_reads(dbg, 0, extension_size, threads, readLogger, false, true, false);
    io::SeqReader reader(reads_lib);
    readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
    coverageStats(logger, dbg);
    if(debug) {
        PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, true);
    }
    dbg.printFastaOld(dir / "final_dbg.fasta");
    printDot(dir / "final_dbg.dot", Component(dbg), readStorage.labeler());
    printGFA(dir / "final_dbg.gfa", Component(dbg), true);
    SaveAllReads(dir/"final_dbg.aln", {&readStorage, &extra_reads});
    readStorage.printReadFasta(logger, dir / "corrected_reads.fasta");
    return {{"corrected_reads", dir/"corrected_reads.fasta"}, {"final_dbg", dir / "final_dbg.gfa"}, {"final_aln", dir / "final_dbg.aln"}};
}

class NoCorrectionStage : public Stage {
public:
    NoCorrectionStage() : Stage(AlgorithmParameters({"k-mer-size=", "window="}, {}, ""),
                                {"reads", "pseudo_reads", "paths"}, {"corrected_reads", "final_dbg", "final_aln"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                      const std::experimental::filesystem::path &dir, bool debug,
                      const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
        size_t w = std::stoi(parameterValues.getValue("window"));
        bool load = false;
        return std::move(NoCorrection(logger, threads, dir, input.find("reads")->second, input.find("pseudo_reads")->second,
                                      input.find("paths")->second, k, w, debug, load));
    }
};
