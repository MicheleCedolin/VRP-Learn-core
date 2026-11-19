add_executable(vrp_tests
  test_main.cpp
  test_methods.cpp)
target_link_libraries(vrp_tests PRIVATE vrplearn)
target_compile_options(vrp_tests PRIVATE ${VRP_FLAGS})

set(REPO_ROOT ${CMAKE_SOURCE_DIR}/..)

# Reduced test set of the tuning repository.  The full project (private) also has the unit tests of the readers, constructions, local
# search and metaheuristics, the frozen reference tables (bit-for-bit comparison) and the Node.js tests of the browser engine; they are not
# part of this repository.  The tests below need data/external for the two benchmark tests (they skip themselves otherwise).
add_test(NAME unit COMMAND vrp_tests WORKING_DIRECTORY ${REPO_ROOT})
