# Optimal control and mathematical modeling of the spread of COVID-19 with multiple delays and among different age groups
The following is an implementation of different optimal control strategies to limit the transmission of the COVID-19 virus, using a SEIHR (Susceptibe, Exposed without symptoms, Infected, Hospitalized, and Recovered) epidemic model. First, we investigated an optimal control strategy which aims at minimizing the number of infected people with and without symptoms and hospitalized individuals (Strategy A). Then, we introduced a more realistic controlled model by incorporating time delay which
specifies the time needed before control actions are initiated (Strategy B). As such, this optimal control strategy includes delay in state and control variables. Finally, we considered another optimal control approach of the model with age groups in control variables whose goal is to describe the dynamics of all age groups population (Strategy C).

This repository contains the data and scripts comprising the three strategies.



<h2>Data and code files</h2>
<h3>A. Optimal control </h3>
The code of this section provides simulations for each compartment (S, E, I, H, R) without control and with the following controls:

1. control u: Sensitization and prevention

2. control v: Quarantine

3. control w: Diagnosis and monitoring

To run this code, go the folder "optimal-control-simple", and run the script "main_optimal_control.m".  

<h3>B. Optimal control with delay</h3>
1. In order to compute optimal control with delay for delay states T=0, T=10 and T=25, go to the folder "optimal-control-delay", and run the script "main_delay.m". After running it, all SEIHR variables and controls u, v, w and z for the three delay states will be saved in a mat file named "covid_retard_SEIHR_and_controls.mat".

2. In order to plot the evolution of each compartment/control for all delay states, run the script "plot_delay.m".

<h3>C. Optimal control with age groups</h3>
In this section, the provided code runs simulations for each compartment (S,E,I,H,R) with and without controls (u and v) with respect to three different age groups:
 
1. Age group 1 : less than 25 years old

2. Age group 2 : between 25 and 60 years old

3. Age group 3 : more than 60 years old

, then plots the evolution of each compartment/control for each and every age group (mentioned above).

For this, go to the folder "optimal-control-agegroups", and run the script "main_agegroups.m".  
