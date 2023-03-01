
 %% 1st step to compute without control and with control with delay T=0 (taux1=0,taux2=0, taux3=0, taux4=0)  
taux1=0;
taux2=0;
taux3=0; 
taux4=0;

[z1,z2,z3,z4,z5,S,E,I,H,R,u11,v1,w1,z11,tt] = compute_optimal_control_delay(taux1,taux2,taux3,taux4);

S_all = [z1;S;zeros(1,length(tt));zeros(1,length(tt))];
E_all = [z2;E;zeros(1,length(tt));zeros(1,length(tt))];
I_all = [z3;I;zeros(1,length(tt));zeros(1,length(tt))];
H_all = [z4;H;zeros(1,length(tt));zeros(1,length(tt))];
R_all = [z5;R;zeros(1,length(tt));zeros(1,length(tt))];
U_all = [u11;zeros(1,length(tt));zeros(1,length(tt))];
V_all = [v1;zeros(1,length(tt));zeros(1,length(tt))];
W_all = [w1;zeros(1,length(tt));zeros(1,length(tt))];        
Z_all = [z11;zeros(1,length(tt));zeros(1,length(tt))];


 %% 2nd step to compute without control and with control with delay T=10 (taux1=7,taux2=8, taux3=9, taux4=10)
taux1=7;
taux2=8;
taux3=9; 
taux4=10;

[z1,z2,z3,z4,z5,S,E,I,H,R,u11,v1,w1,z11,tt] = compute_optimal_control_delay(taux1,taux2,taux3,taux4);

S_all = S_all + [zeros(1,101);zeros(1,101);S(1:101);zeros(1,101)];
E_all = E_all + [zeros(1,101);zeros(1,101);E(1:101);zeros(1,101)];
I_all = I_all + [zeros(1,101);zeros(1,101);I(1:101);zeros(1,101)];
H_all = H_all + [zeros(1,101);zeros(1,101);H(1:101);zeros(1,101)];
R_all = R_all + [zeros(1,101);zeros(1,101);R(1:101);zeros(1,101)];
U_all = U_all + [zeros(1,101);u11(1:101);zeros(1,101)];
V_all = V_all + [zeros(1,101);v1(1:101);zeros(1,101)];
W_all = W_all + [zeros(1,101);w1(1:101);zeros(1,101)];
Z_all = Z_all + [zeros(1,101);z11(1:101);zeros(1,101)];

 %% 2nd step to compute without control and with control with delay  T=25 (taux1=22,taux2=23, taux3=24, taux4=25)
taux1= 22; %12 ./6; % 12 days -> taux1 = 2
taux2= 23; %24 ./6; % 24 days -> taux2 = 4
taux3= 24; % 36./6; % 36 days -> taux3 = 6
taux4= 25; % 48 ./6; % 12 days -> taux4 = 8

[z1,z2,z3,z4,z5,S,E,I,H,R,u11,v1,w1,z11,tt] = compute_optimal_control_delay(taux1,taux2,taux3,taux4);

S_all = S_all + [zeros(1,101);zeros(1,101);zeros(1,101);S(1:101)];
E_all = E_all + [zeros(1,101);zeros(1,101);zeros(1,101);E(1:101)];
I_all = I_all + [zeros(1,101);zeros(1,101);zeros(1,101);I(1:101)];
H_all = H_all + [zeros(1,101);zeros(1,101);zeros(1,101);H(1:101)];
R_all = R_all + [zeros(1,101);zeros(1,101);zeros(1,101);R(1:101)];
U_all = U_all + [zeros(1,101);zeros(1,101);u11(1:101)];
V_all = V_all + [zeros(1,101);zeros(1,101);v1(1:101)];
W_all = W_all + [zeros(1,101);zeros(1,101);w1(1:101)];
Z_all = Z_all + [zeros(1,101);zeros(1,101);z11(1:101)];

save('covid_retard_SEIHR_and_controls.mat','tt','S_all','E_all','I_all','H_all','R_all','U_all','V_all','W_all','Z_all')
