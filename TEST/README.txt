To test the fluidloop program:

1) copy fluidloop.f90 in this folder. In a Bash terminal, get in the loop directory, then:
    cd TEST
    cp ../fluiloop.f90 .
2) compile it. For example, if your fortran compiler is gfortran, execute
    gfortran fluidloop.f90 -o fluidloop.exe
3) execute the fluidloop program
    ./fluidloop.exe
4) check that the output files are identical to reference outputs:
    diff TEST_config.txt TEST_config.txt_ref
    diff TEST_out_0D.txt TEST_out_0D.txt_ref
    diff TEST_out_1D.txt TEST_out_1D.txt_ref
If you did not find any difference, the program works fine.

