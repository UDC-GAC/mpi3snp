#include "Definitions.h"
#include "IOMpi.h"
#include "Search.h"
#include "args.hxx"

Search *configure_search(int &argc, char **&argv) {
    args::ArgumentParser parser(MPI3SNP_DESCRIPTION, MPI3SNP_HELP);
    args::Group required(parser, "Required arguments", args::Group::Validators::All);
    args::Positional<std::string> r_tped(required, "tped-file", "path to TPED file");
    args::Positional<std::string> r_tfam(required, "tfam-file", "path to TFAM file");
    args::Positional<std::string> r_output(required, "output-file", "path to output file");
#ifdef MPI3SNP_USE_GPU
    args::Group gpu_opt(parser, "GPU runtime configuration", args::Group::Validators::DontCare);
    args::ValueFlagList<unsigned int> gpu_ids(gpu_opt, "gpu-ids", "list of GPU ids to use", {'g', "gpu-id"});
#else
    args::Group cpu_opt(parser, "CPU runtime configuration", args::Group::Validators::DontCare);
    args::ValueFlag<unsigned int> cpu_threads(cpu_opt, "thread-num", "number of threads to use per process",
                                              {'t', "threads"});
#endif
    args::Group verb(parser, "Verbosity level", args::Group::Validators::AtMostOne);
    args::Flag verb_b(verb, "benchmarking", "print runtimes", {"benchmarking"});
    args::Flag verb_d(verb, "debug", "print debug information", {"debug"});
    args::ValueFlag<unsigned int> output_num(parser, "output-num", "number of triplets to print in the output file",
                                             {'n', "num-out"});
    args::ValueFlag<bool> use_mi(parser, "mutual-info", "use Mutual Information (the alternative is Information Gain,"
            " default = 1)", {"mi"}, true);
    args::VersionFlag version(parser, "version", "output version information and exit", {'V', "version"});
    args::HelpFlag help(parser, "help", "display this help and exit", {'h', "help"});

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
        IOMpi::Instance().mprint<IOMpi::N>(parser.Help());
        return nullptr;
    }
    catch (args::Version) {
        IOMpi::Instance().mprint<IOMpi::N>(std::string(MPI3SNP_NAME) + " " + MPI3SNP_VERSION + "\n");
        IOMpi::Instance().mprint<IOMpi::N>(std::string(MPI3SNP_LICENSE) + "\n");
        IOMpi::Instance().mprint<IOMpi::N>(std::string("\n") + MPI3SNP_AUTHOR);
        return nullptr;
    }
    catch (const args::ParseError &e) {
        IOMpi::Instance().smprint<IOMpi::N>(std::cerr, e.what());
        IOMpi::Instance().smprint<IOMpi::N>(std::cerr, parser.Help());
        return nullptr;
    }
    catch (const args::ValidationError &e) {
        IOMpi::Instance().smprint<IOMpi::N>(std::cerr, e.what());
        IOMpi::Instance().smprint<IOMpi::N>(std::cerr, parser.Help());
        return nullptr;
    }

    if (verb_d) {
        IOMpi::Instance().Set_print_level(IOMpi::D);
    } else if (verb_b) {
        IOMpi::Instance().Set_print_level(IOMpi::B);
    }

    Search::Builder builder = Search::Builder(args::get(r_tped), args::get(r_tfam), args::get(r_output));

#ifdef MPI3SNP_USE_GPU
    if (gpu_ids) {
        builder.Set_gpu_ids(args::get(gpu_ids));
    }
#else
    if (cpu_threads) {
        builder.Set_cpu_threads(args::get(cpu_threads));
    }
#endif

    if (output_num) {
        builder.Set_num_outputs(args::get(output_num));
    }
    if (use_mi) {
        builder.Set_use_mi(args::get(use_mi));
    }
    return builder.Create_object();
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    double stime, etime;

    /*get the startup time*/
    stime = MPI_Wtime();

    Search *search = configure_search(argc, argv);
    if (search == nullptr) {
        MPI_Finalize();
        return 0;
    }
    search->execute();

    etime = MPI_Wtime();
    IOMpi::Instance().mprint<IOMpi::B>("Overall time: " + std::to_string(etime - stime) + " seconds\n");

    delete search;

    MPI_Finalize();
    return 0;
}
