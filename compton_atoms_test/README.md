This is a test run for the full Monte-Carlo ray-tracing simulation for X-ray photons from AGN. We fire off 6e4 photons, which should take approximately 30 seconds. The time is reported at the end. The structure of the directories is as follows:
* full_test/full_test/ : all the .cpp and .hpp files to run the code, requires external (Boost) libraries
* CMakeLists.txt : setup for compiling via CMake
* build : build directory after running cmake

If you want to test writing to a file, can turn "writeout == true" in full_test/full_test/full_test.cpp.
