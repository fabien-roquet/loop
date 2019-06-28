** Instructions to test the fluidloop program (version 1)

1. copy fluidloop.f90 in your working directory. In a Bash terminal, get in the loop directory, then:

    ```
    WORKDIR=my_working_directory
    cp fluiloop.f90 $WORKDIR
    cd $WORKDIR
    ```
    
2. compile it. For example, if your fortran compiler is gfortran, execute

    ```
    gfortran fluidloop.f90 -o fluidloop.exe
    ```
    
3. execute the fluidloop program

    ```
    ./fluidloop.exe
    ```
    
    Important: If output files already exists in your working directory, the program will stop. Rename or remove them before running the loop program.
    
4. check that the output files are identical to reference outputs:

    ```
    diff TEST_config.txt TEST_config.txt_ref
    diff TEST_out_0D.txt TEST_out_0D.txt_ref
    diff TEST_out_1D.txt TEST_out_1D.txt_ref
    ```
    
If you did not find any difference, it means that the program works fine.

