function plot_group_age(y,z1,z2,z3,z4,z5,S,E,I,H,R,u11,v1)

%% plot the evolution of the population size without control
figure;
plot(y(1,:),z1(1,:),'r',y(1,:),z1(2,:),'g',y(1,:),z1(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Susceptible');
legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years');
hold off
grid
title('Population size without control');

%% plot the evolution of the population size with control
figure;
plot(y(1,:),S(1,:),'r',y(1,:),S(2,:),'g',y(1,:),S(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Susceptible');
legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years');
hold off
grid
title('Population size with control');


%% plot the evolution of the population E without control
figure;
plot(y(1,:),z2(1,:),'r',y(1,:),z2(2,:),'g',y(1,:),z2(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Exposed');
h=legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years','Location','northeast');
set(h,'Interpreter','latex')
grid on;
title('Exposed without control (a)');
maximumvalue = max(z2(:));

%% plot the evolution of the population E with control
figure;
plot(y(1,:),E(1,:),'r',y(1,:),E(2,:),'g',y(1,:),E(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Exposed');
h=legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years','Location','northeast');
set(h,'Interpreter','latex')
grid on;

title('Exposed with control (b)');
ylim([0 maximumvalue+1])

%% plot the evolution of the population I without control
figure;
plot(y(1,:),z3(1,:),'r',y(1,:),z3(2,:),'g',y(1,:),z3(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Infected');
h=legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years','Location','northeast');
set(h,'Interpreter','latex')
grid on;
title('Infected without control (a)');
maximumvalue = max(z3(:));

%% plot the evolution of the population I with control
figure;
plot(y(1,:),I(1,:),'r',y(1,:),I(2,:),'g',y(1,:),I(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Infected');
h=legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years','Location','northeast');
set(h,'Interpreter','latex')
grid on;
title('Infected with control (b)');
ylim([0 maximumvalue+1])

%% plot the evolution of the population H without control
figure;
plot(y(1,:),z4(1,:),'r',y(1,:),z4(2,:),'g',y(1,:),z4(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Hospitalized');
h=legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years','Location','northeast');
set(h,'Interpreter','latex')
grid on;
title('Hospitalized  without control (a)');
maximumvalue = max(z4(:));

%% plot the evolution of the population H with control
figure;
plot(y(1,:),H(1,:),'r',y(1,:),H(2,:),'g',y(1,:),H(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Hospitalized');
h=legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years','Location','northeast');
set(h,'Interpreter','latex')
grid on;

title('Hospitalized with control (b)');
ylim([0 maximumvalue+1])

%% plot the evolution of the population R without control
figure;
plot(y(1,:),z5(1,:),'r',y(1,:),z5(2,:),'g',y(1,:),z5(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Recovred');
h=legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years','Location','northeast');
set(h,'Interpreter','latex')
grid on;
title('Recovered without control (a)');
maximumvalue = max(z5(:));

%% plot the evolution of the population R with control
figure;
plot(y(1,:),R(1,:),'r',y(1,:),R(2,:),'g',y(1,:),R(3,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Recovred');
h=legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years','Location','northeast');
set(h,'Interpreter','latex')
grid on;

title('Recovered with control (b)');
ylim([0 maximumvalue+1])

%% plot the evolution of the control u
figure;
hold on
plot(y(1,:),u11(1,:),'r',y(1,:),u11(2,:),'g',y(1,:),u11(3,:),'b');
legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years');
xlabel('time (in days)');
ylabel('Control u');
title('(u)');
%axis([0 15 0 1]);
grid on;

%% plot the evolution of the control v
figure;
hold on
plot(y(1,:),v1(1,:),'r',y(1,:),v1(2,:),'g',y(1,:),v1(3,:),'b');
legend('Age group under 25 years','Age group between 25 and 65 years', 'Age group over the 65 years');
xlabel('time (in days)');
ylabel('Control v');
title('(v)');
%axis([0 15 0 1]);
grid on;

end