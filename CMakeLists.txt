#add_executable(dc.darts.exe adc.output.darts.cpp dc.output.darts.cpp extbuild.output.darts.cpp rbt.output.darts.cpp jobcntl.output.darts.cpp c_print_results.output.darts.cpp c_timers.output.darts.cpp wtime.output.darts.cpp ../../darts-lib/utils.cpp)
add_executable(lu.darts.exe lu.output.darts.cpp c_print_results.output.darts.cpp c_timers.output.darts.cpp wtime.output.darts.cpp ../../darts-lib/utils.cpp)
target_link_libraries(lu.darts.exe darts)
find_package(OpenMP)

#add_executable(lu.exe adc.c dc.c extbuild.c rbt.c jobcntl.c c_print_results.c c_timers.c wtime.c)
LINK_LIBRARIES(m)
add_executable(lu.exe lu.c c_print_results.c c_timers.c wtime.c)
target_compile_definitions(lu.exe PRIVATE "-D_OPENMP")
target_compile_options(lu.exe PRIVATE "-fopenmp")
target_link_libraries(lu.exe PRIVATE OpenMP::OpenMP_C)

include_directories(../../darts-lib)
