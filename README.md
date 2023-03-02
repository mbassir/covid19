# covid19
Optimal control and mathematical modeling of the spread of COVID-19 with multiple delays and among different age groups

<h2>Data and code files</h2>
<h3>A. Optimal control with age groups </h3>

<h3>B. Optimal control with age groups </h3>
In this section, the provided code runs simulations for each compartment (S,E,I,H,R) with and without controls (u and v) with respect to three different age groups:

1. Age group 1 : less than 25 years old

2. Age group 2 : between 25 and 60 years old

3. Age group 3 : more than 60 years old

, then plots the evolution of each compartment/control for each and every age group (mentioned above).

For this, go to the folder "optimal-control-agegroups", and run the script "main_agegroups.m".  

<h3>C. Optimal control with delay</h3>
1. In order to compute optimal control with delay for delay states T=0, T=10 and T=25, go to the folder "optimal-control-delay", and run the script "main_delay.m". After running it, all SEIHR variables and controls u, v, w and z for the three delay states will be saved in a mat file named "covid_retard_SEIHR_and_controls.mat".

2. In order to plot the evolution of each compartment/control for all delay states, run the script "plot_delay.m".
