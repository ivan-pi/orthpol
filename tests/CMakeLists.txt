project(orthpol)

link_directories(${CMAKE_BINARY_DIR}/lib)

macro(ADDTEST name)
    add_executable(${name} ${name}.f)
    target_link_libraries(${name} p::orthpol)
    add_test(${name} ${PROJECT_BINARY_DIR}/${name})
endmacro(ADDTEST)

macro(ADDTESTF90 name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} p::orthpol)
    add_test(${name} ${PROJECT_BINARY_DIR}/${name})
endmacro(ADDTESTF90)

ADDTESTF90(test1)

#ADDTEST(test1)
ADDTEST(test2)
ADDTEST(test3)
ADDTEST(test4)
ADDTEST(test5)
ADDTEST(test6)
ADDTEST(test7)
ADDTEST(test8)
ADDTEST(test9)
ADDTEST(test10)
ADDTEST(test11)