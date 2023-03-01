load('covid_retard_SEIHR_and_controls.mat');

tt= tt(1:101);

%% plot the evolution of the population S
figure;
plot(tt,S_all(1,:),'r',tt,S_all(2,:),'g',tt,S_all(3,:),'b',tt,S_all(4,:),'y');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Susceptible');
legend('S without control','S with control and without delay', 'S with control and with delay T=10', 'S with control and with delay T=25');
hold off
grid

%% plot the evolution of the population E
figure;
plot(tt,E_all(1,:),'r',tt,E_all(2,:),'g',tt,E_all(3,:),'b',tt,E_all(4,:),'y');
hold on
title('(E)');
xlabel('Time(days)');  ylabel('Cumulative number of Exposed');
legend('E without control','E with control and without delay', 'E with control and with delay T=10', 'E with control and with delay T=25');
hold off
grid


%% plot the evolution of the population E
figure;
plot(tt,I_all(1,:),'r',tt,I_all(2,:),'g',tt,I_all(3,:),'b',tt,I_all(4,:),'y');
hold on
title('(I)');
xlabel('Time(days)');  ylabel('Cumulative number of Infected');
legend('I without control','I with control and without delay', 'I with control and with delay T=10', 'I with control and with delay T=25');
hold off
grid

%% plot the evolution of the population H
figure;
plot(tt,H_all(1,:),'r',tt,H_all(2,:),'g',tt,H_all(3,:),'b',tt,H_all(4,:),'y');
hold on
title('(H)');
xlabel('Time(days)');  ylabel('Cumulative number of Hospitalized');
legend('H without control','H with control u and without delay', 'H with control u and with delay T=10', 'H with control u and with delay T=25');
hold off
grid


%% plot the evolution of the population R
figure;
plot(tt,R_all(1,:),'r',tt,R_all(2,:),'g',tt,R_all(3,:),'b',tt,R_all(4,:),'y');
hold on
title('(R)');
xlabel('Time(days)');  ylabel('Cumulative number of Recovred');
legend('R without control','R with control and without delay', 'R with control and with delay T=10', 'R with control and with delay T=25');
hold off
grid

%% plot the evolution of the control u
figure;
plot(tt,U_all(1,:),'r',tt,U_all(2,:),'g',tt,U_all(3,:),'b');
hold on
title('(u)');
xlabel('Time(days)');  ylabel('Control u');
legend('u without delay', 'u with delay T=10','u with delay T=25');
hold off
grid


%% plot the evolution of the control v
figure;
plot(tt,V_all(1,:),'r',tt,V_all(2,:),'g',tt,V_all(3,:),'b');
hold on
title('(v)');
xlabel('Time(days)');  ylabel('Control v');
legend('v without delay', 'v with delay T=10','v with delay T=25');
hold off
grid

%% plot the evolution of the control w
figure;
plot(tt,W_all(1,:),'r',tt,W_all(2,:),'g',tt,W_all(3,:),'b');
hold on
title('(w)');
xlabel('Time(days)');  ylabel('Control w');
legend('w without delay', 'w with delay T=10','w with delay T=25');
hold off
grid


%% plot the evolution of the control z
figure;
plot(tt,Z_all(1,:),'r',tt,Z_all(2,:),'g',tt,Z_all(3,:),'b');
hold on
title('(z)');
xlabel('Time(days)');  ylabel('Control z');
legend('z without delay', 'z with delay T=10','z with delay T=25');
hold off
grid



% end
