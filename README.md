# KMCcGaN
Kinetic Monte Carlo simulations of c-GaN growth
* This code KMC code is used to simulate epitaxial growth on GaN(0001) surface for both flat and vicinal surface.*
The list of files and their descriptions is:
1. main.c - Main program
2. variables.h - Values of parameters and the declaration of global variables of the program
3. header.h - List of subroutines used
4. Subroutine files: 
AtomSelection.c
countNeighb.c
External_Flux.c
H3_AdGa.c
H3_N.c
H3_N_rot.c
H3_sbN.c
H4_Ga.c
H7_AdGa.c
H7_N.c
H7_N_rot.c
H7_sbN.c
H8_Ga.c
Height_Calculation.c
Initial_Conditions.c
move_AdGa.c
move_Ga.c
move_N.c
move_sbN.c
OutGaNfiles.c
ReadInput.c
ReadOutput.c
stable_AdGa.c
stable_Ga.c
stable_N.c
stable_N_OGa.c



The steps to execute the code are : 

1. Edit Input variables: All variables are in the file variables.h. Here it is possible to decide the length of simulation, the relative flux ratio of Ga to N, the overall flux and the temperature.  All the allowed process and their energy barriers are listed in variables.h. 
 In case the simulation is restarted from some previous simulation, the value of outNfiles in the main.c should be edited.
 In case this code is used for vicinal surface: an if..else statement has to be uncommmented in move_Ga.c, move_N.c, move_AdGa.c, move_sbN.c, countNeigh.c 
 
2. Edit compile.sh  so the correct export path of gsl libraries is provided. Run ./compile.sh

3. Execute the code using ./a.out After the code is executed: It generates one file named input and several files named OutputGaN%t which contain output files.

4. These OutputGaN files are used for post processing such as calculating roughness and making snapshots.


This code was written by Manjusha Chugh, Razia and Madhav Ranganathan at IIT Kanpur. The current version is December 12, 2024. Any comments/suggestions can be commuicated to madhavr@iitk.ac.in 
