# ECM3446 Calculating and Parallelising a Numerical Solution to a 2D Advection Equation

## Introduction
The intention of this project is to calculate a numerical solution to a two dimensional advection equation that simulates the movement of clouds of material in the atmospheric boundary layer. Furthermore, this program uses OpenMP to parallelise parts of the program to increase efficiency.

## Prerequisites and Installation

This project uses the C programming languages as its base with GNUplot for plotting the graphs. Additionally the paralleisation extension OpenMP is used to parallelise the program and speed up its execution. Both GNUplot and OpenMP are included in most systems (if your system uses an up to date GCC compiler and has Python installed you are in the clear). Otherwise consult the following pages for installation instructions:

    - http://www.gnuplot.info/download.html

    - https://www.openmp.org/resources/openmp-compilers-tools/

## Project Tutorial
The source code for the final version of the program is stored in the src directory, including the initial.dat and final.dat data files used to plot the graphs stored in the plots directory. Additionally, the gnuplot program used to plot the vertically averaged distribution of the advection program and the plot_vertically_averaged_distribution.c file used to plot this graph are included.

Only the final version of the program is provided so recreation of the graphs required by question 2.2 of the spec can only come from reverting the changes made in 2.3 (specification can be found in the docs directory), recompiling and rexecuting the program.

For the steps below you must navigate to the src directory of the project folder (cd ../src).

The main advection2D.c program can be compiled and executed as follows:

    - gcc advection2D.c -o advection2D -fopenmp -std=c99 -lm
    - ./advection2D
    
To create a non-parallelised version of this program just compile without the OpemMP compiler flag as such:

    - gcc advection2D.c -o advection2D -std=c99 -lm
    - ./advection2D
    
A separate gnuplot file was used to plot the graph for the vertically averaged distribution, this file can be run via:

    - gnuplot "plot_vertically_averaged_distribution.c"
    
The graphs for question 2.3 and 2.2 (after reverting the changes required by 2.3 and recompiling and re-executing the program) can be created via the following command:

    - gnuplot "plot.c"

## Testing
No testing was required in the spec, and so no testing was implemented.

### Notes
All functionality that was required by the specification was provided.

## Details

#### Authors
Benjamin Theron

#### License
MIT License

