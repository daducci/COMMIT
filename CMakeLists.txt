cmake_minimum_required(VERSION 3.21...4.0)
project(dmri-commit CXX C)
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
include(UseCython)

set(CYTHON_ARGS "-3")

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
        "import sys; sys.path.insert(0,'${CMAKE_CURRENT_SOURCE_DIR}'); from setup_operator import write_operator_c_file; write_operator_c_file()"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMAND_ERROR_IS_FATAL ANY
)

# trx-cpp via FetchContent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
include(FetchContent)

set(TRX_BUILD_TESTS      OFF)
set(TRX_BUILD_EXAMPLES   OFF)
set(TRX_BUILD_DOCS       OFF)
set(TRX_BUILD_BENCHMARKS OFF)
set(TRX_ENABLE_NIFTI     OFF)
set(TRX_ENABLE_INSTALL   OFF)

FetchContent_Declare(
    trx-cpp
    GIT_REPOSITORY https://github.com/tee-ar-ex/trx-cpp.git
    GIT_TAG        main
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(trx-cpp)

# json11.cpp is missing '#include <cstdint>'; force-include it
if(MSVC)
    target_compile_options(trx PRIVATE /FIcstdint)
else()
    target_compile_options(trx PRIVATE -include cstdint)
endif()

cython_transpile(commit/trk2dictionary/trk2dictionary.pyx LANGUAGE CXX OUTPUT_VARIABLE trk2dict_src)
Python_add_library(trk2dictionary MODULE "${trk2dict_src}" WITH_SOABI)
target_compile_features(trk2dictionary PRIVATE cxx_std_17)
target_compile_options(trk2dictionary PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-w>)
target_include_directories(trk2dictionary PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/commit/trk2dictionary")
target_link_libraries(trk2dictionary PRIVATE trx-cpp::trx Python::NumPy)
install(TARGETS trk2dictionary DESTINATION commit)

cython_transpile(commit/core.pyx LANGUAGE CXX OUTPUT_VARIABLE core_src)
Python_add_library(core MODULE "${core_src}" WITH_SOABI)
target_compile_features(core PRIVATE cxx_std_11)
target_compile_options(core PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-w>)
target_link_libraries(core PRIVATE Python::NumPy)
install(TARGETS core DESTINATION commit)

cython_transpile(commit/proximals.pyx LANGUAGE CXX OUTPUT_VARIABLE proximals_src)
Python_add_library(proximals MODULE "${proximals_src}" WITH_SOABI)
target_compile_features(proximals PRIVATE cxx_std_11)
target_compile_options(proximals PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-w>)
target_link_libraries(proximals PRIVATE Python::NumPy)
install(TARGETS proximals DESTINATION commit)

cython_transpile(commit/models.pyx LANGUAGE CXX OUTPUT_VARIABLE models_src)
Python_add_library(models MODULE "${models_src}" WITH_SOABI)
target_compile_features(models PRIVATE cxx_std_11)
target_compile_options(models PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-w>)
target_link_libraries(models PRIVATE Python::NumPy)
install(TARGETS models DESTINATION commit)

cython_transpile(commit/operator/operator.pyx LANGUAGE C OUTPUT_VARIABLE operator_src)
Python_add_library(commit_operator MODULE
    "${operator_src}"
    "${CMAKE_CURRENT_SOURCE_DIR}/commit/operator/operator_c.c"
    WITH_SOABI
)
set_target_properties(commit_operator PROPERTIES OUTPUT_NAME operator)
target_compile_options(commit_operator PRIVATE
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-w -O3 -Ofast>
    $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /DHAVE_STRUCT_TIMESPEC>
)
target_link_libraries(commit_operator PRIVATE Python::NumPy)
if(WIN32)
    target_include_directories(commit_operator PRIVATE "$ENV{PTHREAD_WIN_INCLUDE}")
    target_link_directories(commit_operator PRIVATE "$ENV{PTHREAD_WIN_LIB}")
    target_link_libraries(commit_operator PRIVATE pthread)
endif()
install(TARGETS commit_operator DESTINATION commit/operator)
