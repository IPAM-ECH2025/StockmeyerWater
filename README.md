## StockmeyerWater

This repositoty contains applications directories for the PRISMS-PF library. 

In order to run, in your computer of HPC cluster, you must install de [deal.II](https://dealii.org/current_release/download/) library as well as the [PRISMS-PF](https://prisms-center.github.io/phaseField/doxygen/3.0.0-pre/installation.html) library. Alternatively, you can install the PRISMS-PF [Docker](https://prisms-center.github.io/phaseField/doxygen/3.0.0-pre/docker.html) version, which already contains deal.II. 

Once these two libraries are installed, you will need to set the environment variable that specifies where the PRISMS-PF library is installed: 
```bash
export PRISMS_PF_DIR=</path/to/installation/>
```
After that, you should be able to compile each application by typin the following within the application folder.
```bash
cmake . && make
```
This will generate two executable files: **main-release** and **main-debug**. Debug and release are compiler configurations. Debug mode is slower, but contains fewer optimiziations and more meaningful error messages. This makes it ideal for application/model code development. Release mode has fewer "safety features" and meaningful error messages but more optimizations (faster runtime).

Debug execution (serial runs):
```bash
$ ./main-debug
```
Release execution (parallel runs):
```bash
$ mpirun -np <nprocs> ./main-release
```
Here, `<nprocs>` denotes the number of parallel tasks you want to use run the code.