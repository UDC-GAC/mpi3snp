#include "Definitions.h"
#include "IOMpi.h"
#include "Search.h"
#include "args.hxx"

Search *configure_search(int &argc, char **argv){
    args::ArgumentParser parser(MPI3SNP_DESCRIPTION, MPI3SNP_HELP);
    args::Group group(parser, "Required arguments", args::Group::Validators::All);
    args::Positional<std::string> tped(group, "tped-file", "path to TPED file");
    args::Positional<std::string> tfam(group, "tfam-file", "path to TFAM file");
    args::Positional<std::string> output(group, "output-file", "path to output file");
    args::ValueFlagList<unsigned int> gpu_ids(parser, "gpu-ids", "list of GPU ids to use", {'g', "gpu-id"});
    args::ValueFlag<unsigned int> output_num(parser, "output-num", "number of triplets to print in the output file",
                                             {'n', "num-out"});
    args::ValueFlag<bool> use_mi(parser, "mutual-info", "use Mutual Information (the alternative is Information Gain,"
            " default = 1)", {"mi"}, true);
    args::Group verb(parser, "Verbosity level", args::Group::Validators::AtMostOne);
    args::Flag verb_b(verb, "benchmarking", "print runtimes", {"benchmarking"});
    args::Flag verb_d(verb, "debug", "print debug information", {"debug"});
    args::VersionFlag version(parser, "version", "output version information and exit", {'V', "version"});
    args::HelpFlag help(parser, "help", "display this help and exit", {'h', "help"});

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
        std::cout << parser;
        return nullptr;
    }
    catch (args::Version) {
        std::cout << MPI3SNP_NAME << " " << MPI3SNP_VERSION << std::endl;
        std::cout << MPI3SNP_LICENSE << std::endl;
        std::cout << std::endl << MPI3SNP_AUTHOR << std::endl;
        return nullptr;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return nullptr;
    }
    catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return nullptr;
    }

    if (verb_b){
        IOMpi::Instance().Set_print_level(IOMpi::B);
    }
    if (verb_d){
        IOMpi::Instance().Set_print_level(IOMpi::D);
    }

    Search::Builder builder = Search::Builder(args::get(tped), args::get(tfam), args::get(output));
    if (gpu_ids) {
        builder.Set_gpu_ids(args::get(gpu_ids));
    }
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

    // TODO: 1.Select appropiate CPU/GPU implementation 2.Adequate arguments
    Search *search = configure_search(argc, argv);
    if (search == nullptr){
        IOMpi::Instance().Deallocate_MPI_resources();
        MPI_Finalize();
        return 0;
    }
    search->execute();

    etime = MPI_Wtime();
    IOMpi::Instance().Mprintf<IOMpi::B>("Overall time: %.2f seconds\n", etime - stime);

    delete search;

    IOMpi::Instance().Deallocate_MPI_resources();
    MPI_Finalize();
    return 0;
}
