project(ARGS)
set(ARGS_SOURCE_FILES args.hxx)
add_library(ARGS INTERFACE)
target_sources(ARGS INTERFACE args.hxx)