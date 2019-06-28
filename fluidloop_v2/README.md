** Instructions to test the fluidloop program (version 2)

1. You need to install a fortran compiler (e.g. gfortran). You need also to install 
netcdf and netcdf-fortran with the same compiler.

2. Copy the source files *.f90 and the Makefile in your working directory.
Modify the Makefile to point to the right fortran compiler and netcdf libraries.

    ```
    cp *.f90 Makefile $WORKDIR
    cd $WORKDIR
    ```


2. Compile the code with the command make.

    ```
    make
    make clean
    ```
    
3. execute the fluidloop program

    ```
    ./fluidloop
    ```
    
    Important: If output files already exists in your working directory, the program will stop. Rename or remove them before running the loop program.
    
4. check that the output files are identical to reference outputs:

    ```
    diff TEST_config.txt TEST_config.txt_ref
    diff TEST_output.nc TEST_output.nc_ref
    ```
    
If you did not find any difference, it means that the program works fine.
Note that some difference may arise on the *_output.nc files depending on the version of 
the netcdf libraries you are using. A safer test will be to compare the data in the netcdf file
using the ncdiff command (nco toolbox) or using a scientific program such as Python or Matlab.

