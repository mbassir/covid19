%%%%%%%%%%%% Optimal control %%%%%%%%%%%%%
% This code provides simulations for each compartment (S, E, I, H, R)
% without control and with three controls:
%               control u: Sensitization and prevention
%               control v: Quarantine
%               control w: Diagnosis and monitoring

clc;
clear all;
global eta beta mu N alpha theta delta1 lamda gamma delta2

eta=10; mu=0.0045; delta1=0.001; delta2=0.004; beta=0.045; alpha=0.022;
theta=0.65; lamda=0.024; gamma=0.45;
A=0.000010; B=0.000010; C=0.0000010;


K=1e+6;

delta=0.001;
S0=8000; I0=170; E0=300; H0=30; R0=20;
S=S0; E=E0; I=I0; H=H0; R=R0;
MaxTime=100;
test=-1;
NN=100;
N=S+E+I+H+R;
t=linspace(0,MaxTime,NN+1);
h=MaxTime/NN;
h2=h/2;

u=zeros(1,NN+1);
v=zeros(1,NN+1);
w=zeros(1,NN+1);

S=zeros(1,NN+1); E=zeros(1,NN+1); I=zeros(1,NN+1); H=zeros(1,NN+1);
R=zeros(1,NN+1);
S(1)=S0; E(1)=E0; I(1)=I0; H(1)=H0; R(1)=R0;

z1=zeros(1,NN+1); z2=zeros(1,NN+1); z3=zeros(1,NN+1); z4=zeros(1,NN+1);
z5=zeros(1,NN+1);
z1(1)=S0; z2(1)=E0; z3(1)=I0; z4(1)=H0; z5(1)=R0;

sigma1=zeros(1,NN+1); sigma2=zeros(1,NN+1); sigma3=zeros(1,NN+1);
sigma4=zeros(1,NN+1); sigma5=zeros(1,NN+1);
%condition de transversalité
sigma2(NN+1)=1; sigma3(NN+1)=1;

x=[S E I H R]'; z=[z1 z2 z3 z4 z5]';

sigma=[sigma1 sigma2 sigma3 sigma4 sigma5]';
u1=u';

for i=1:NN
    
    k11=eta - beta*z1(i)*z2(i)/N - mu*z1(i) ;
    k12=beta*z1(i)*z2(i)/N - (mu+alpha+theta)*z2(i) ;
    k13=alpha*z2(i) - (mu+lamda+delta1)*z3(i) ;
    k14=lamda*z3(i) - (mu+gamma+delta2)*z4(i) ;
    k15=gamma*z4(i) + theta*z2(i) - mu*z5(i) ;
    
    k21=eta - beta*(z1(i)+h2*k11)*(z2(i)+h2*k12)/N - mu*(z1(i)+h2*k11) ;
    k22=beta*(z1(i)+h2*k11)*(z2(i)+h2*k12)/N - (mu+alpha+theta)*(z2(i)+h2*k12) ;
    k23=alpha*(z2(i)+h2*k12) - (mu+lamda+delta1)*(z3(i)+h2*k13) ;
    k24=lamda*(z3(i)+h2*k13) - (mu+gamma+delta2)*(z4(i)+h2*k14) ;
    k25=gamma*(z4(i)+h2*k14) + theta*(z2(i)+h2*k12) - mu*(z5(i)+h2*k15) ;
    
    k31=eta - beta*(z1(i)+h2*k21)*(z2(i)+h2*k22)/N - mu*(z1(i)+h2*k21) ;
    k32=beta*(z1(i)+h2*k21)*(z2(i)+h2*k22)/N - (mu+alpha+theta)*(z2(i)+h2*k22) ;
    k33=alpha*(z2(i)+h2*k22) - (mu+lamda+delta1)*(z3(i)+h2*k23) ;
    k34=lamda*(z3(i)+h2*k23) - (mu+gamma+delta2)*(z4(i)+h2*k24) ;
    k35=gamma*(z4(i)+h2*k24) + theta*(z2(i)+h2*k22) - mu*(z5(i)+h2*k25) ;
    
    k41=eta - beta*(z1(i)+h2*k31)*(z2(i)+h2*k32)/N - mu*(z1(i)+h2*k31) ;
    k42=beta*(z1(i)+h2*k31)*(z2(i)+h2*k32)/N - (mu+alpha+theta)*(z2(i)+h2*k32) ;
    k43=alpha*(z2(i)+h2*k32) - (mu+lamda+delta1)*(z3(i)+h2*k33) ;
    k44=lamda*(z3(i)+h2*k33) - (mu+gamma+delta2)*(z4(i)+h2*k34) ;
    k45=gamma*(z4(i)+h2*k34) + theta*(z2(i)+h2*k32) - mu*(z5(i)+h2*k35) ;
    
    z1(i+1)=z1(i)+(h/6)*(k11+2*k21+2*k31+k41);
    z2(i+1)=z2(i)+(h/6)*(k12+2*k22+2*k32+k42);
    z3(i+1)=z3(i)+(h/6)*(k13+2*k23+2*k33+k43);
    z4(i+1)=z4(i)+(h/6)*(k14+2*k24+2*k34+k44);
    z5(i+1)=z5(i)+(h/6)*(k15+2*k25+2*k35+k45);
    
end

count = 1;

%    while(test<0)

while(test<0 && count<=20)
    
    oldu=u;
    oldu1=u1; oldv=v; oldw=w;
    oldx=x;
    oldsigma=sigma;
    
    for i=1:NN
        
        
        k11=eta - beta*(1-u(i))*S(i)*E(i)/N - mu*S(i) ;
        k12=beta*(1-u(i))*S(i)*E(i)/N - (mu+alpha+theta)*E(i) - w(i)*E(i) ;
        k13=alpha*E(i) - (mu+lamda+delta1)*I(i) - v(i)*I(i) ;
        k14=lamda*I(i) - (mu+gamma+delta2)*H(i) + v(i)*I(i) + w(i)*E(i) ;
        k15=gamma*H(i) + theta*E(i) - mu*R(i) ;
        
        k21=eta - beta*(1-0.5*(u(i)+u(i+1)))*(S(i)+h2*k11)*(E(i)+h2*k12)/N - mu*(S(i)+h2*k11) ;
        k22=beta*(1-0.5*(u(i)+u(i+1)))*(S(i)+h2*k11)*(E(i)+h2*k12)/N - (mu+alpha+theta)*(E(i)+h2*k12) - 0.5*(w(i)+w(i+1))*(E(i)+h2*k12) ;
        k23=alpha*(E(i)+h2*k12) - (mu+lamda+delta1)*(I(i)+h2*k13) - 0.5*(v(i)+v(i+1))*(I(i)+h2*k13) ;
        k24=lamda*(I(i)+h2*k13) - (mu+gamma+delta2)*(H(i)+h2*k14) + 0.5*(v(i)+v(i+1))*I(i) + 0.5*(w(i)+w(i+1))*(E(i)+h2*k12) ;
        k25=gamma*(H(i)+h2*k14) + theta*(E(i)+h2*k12) - mu*(R(i)+h2*k15) ;
        
        k31=eta - beta*(1-0.5*(u(i)+u(i+1)))*(S(i)+h2*k21)*(E(i)+h2*k22)/N - mu*(S(i)+h2*k21) ;
        k32=beta*(1-0.5*(u(i)+u(i+1)))*(S(i)+h2*k21)*(E(i)+h2*k22)/N - (mu+alpha+theta)*(E(i)+h2*k22) - 0.5*(w(i)+w(i+1))*E(i) ;
        k33=alpha*(E(i)+h2*k22) - (mu+lamda+delta1)*(I(i)+h2*k23) - 0.5*(v(i)+v(i+1))*I(i) ;
        k34=lamda*(I(i)+h2*k23) - (mu+gamma+delta2)*(H(i)+h2*k24) + 0.5*(v(i)+v(i+1))*I(i) + 0.5*(w(i)+w(i+1))*E(i) ;
        k35=gamma*(H(i)+h2*k24) + theta*(E(i)+h2*k22) - mu*(R(i)+h2*k25) ;
        
        k41=eta - beta*(1-u(i+1))*(S(i)+h2*k31)*(E(i)+h2*k32)/N - mu*(S(i)+h2*k31) ;
        k42=beta*(1-u(i+1))*(S(i)+h2*k31)*(E(i)+h2*k32)/N - (mu+alpha+theta)*(E(i)+h2*k32) - w(i+1)*(E(i)+h2*k32) ;
        k43=alpha*(E(i)+h2*k32) - (mu+lamda+delta1)*(I(i)+h2*k33) - v(i+1)*(I(i)+h2*k33) ;
        k44=lamda*(I(i)+h2*k33) - (mu+gamma+delta2)*(H(i)+h2*k34) + v(i+1)*I(i) + w(i+1)*(E(i)+h2*k32) ;
        k45=gamma*(H(i)+h2*k34) + theta*(E(i)+h2*k32) - mu*(R(i)+h2*k35) ;
        
        S(i+1)=S(i)+(h/6)*(k11+2*k21+2*k31+k41);
        E(i+1)=E(i)+(h/6)*(k12+2*k22+2*k32+k42);
        I(i+1)=I(i)+(h/6)*(k13+2*k23+2*k33+k43);
        H(i+1)=H(i)+(h/6)*(k14+2*k24+2*k34+k44);
        R(i+1)=R(i)+(h/6)*(k15+2*k25+2*k35+k45);
        
    end
    
    fprintf('counting iterations %f \n',count);
    
    for i=1:NN
        j=NN+2-i;
        
        
        k11=(beta*(1-u(j))*E(j)/N)*(sigma1(j)-sigma2(j)) + mu*sigma1(j);
        k12=-1+(beta*(1-u(j))*S(j)/N)*(sigma1(j)-sigma2(j))+(mu+alpha+theta+w(j))*sigma2(j) - alpha*sigma3(j) - w(j)*sigma4(j) - theta*sigma5(j);
        k13=-1+(mu+alpha+delta1+v(j))*sigma3(j) - (lamda+v(j))*sigma4(j);
        k14=(mu+gamma+delta2)*sigma4(j) - gamma*sigma5(j);
        k15=mu*sigma5(j);
        
        k21=(beta*(1-0.5*(u(j)+u(j-1)))*0.5*(E(j)+E(j-1))/N)*((sigma1(j)-h2*k11)-(sigma2(j)-h2*k12)) + mu*(sigma1(j)-h2*k11);
        k22=-1+(beta*(1-0.5*(u(j)+u(j-1)))*0.5*(S(j)+S(j-1))/N)*((sigma1(j)-h2*k11)-(sigma2(j)-h2*k12))+(mu+alpha+theta+0.5*(w(j)+w(j-1)))*(sigma2(j)-h2*k12) - ...
            - alpha*(sigma3(j)-h2*k13) - 0.5*(w(j)+w(j-1))*(sigma4(j)-h2*k14) - theta*(sigma5(j)-h2*k15);
        k23=-1+(mu+alpha+delta1+0.5*(v(j)+v(j-1)))*(sigma3(j)-h2*k13) - (lamda+0.5*(v(j)+v(j-1)))*(sigma4(j)-h2*k14);
        k24=(mu+gamma+delta2)*(sigma4(j)-h2*k14) - gamma*(sigma5(j)-h2*k15);
        k25=mu*(sigma5(j)-h2*k15);
        
        k31=(beta*(1-0.5*(u(j)+u(j-1)))*0.5*(E(j)+E(j-1))/N)*((sigma1(j)-h2*k21)-(sigma2(j)-h2*k22)) + mu*(sigma1(j)-h2*k21);
        k32=-1+(beta*(1-0.5*(u(j)+u(j-1)))*0.5*(S(j)+S(j-1))/N)*((sigma1(j)-h2*k21)-(sigma2(j)-h2*k22))+(mu+alpha+theta+0.5*(w(j)+w(j-1)))*(sigma2(j)-h2*k22) - ...
            - alpha*(sigma3(j)-h2*k23) - 0.5*(w(j)+w(j-1))*(sigma4(j)-h2*k24) - theta*(sigma5(j)-h2*k25);
        k33=-1+(mu+alpha+delta1+0.5*(v(j)+v(j-1)))*(sigma3(j)-h2*k23) - (lamda+0.5*(v(j)+v(j-1)))*(sigma4(j)-h2*k24);
        k34=(mu+gamma+delta2)*(sigma4(j)-h2*k24) - gamma*(sigma5(j)-h2*k25);
        k35=mu*(sigma5(j)-h2*k25);
        
        k41=(beta*(1-u(j-1))*E(j-1)/N)*((sigma1(j)-h2*k31)-(sigma2(j)-h2*k32)) + mu*(sigma1(j)-h2*k31);
        k42=-1+(beta*(1-u(j-1))*S(j-1)/N)*((sigma1(j)-h2*k31)-(sigma2(j)-h2*k32))+(mu+alpha+theta+w(j-1))*(sigma2(j)-h2*k32) - ...
            - alpha*(sigma3(j)-h2*k33) - w(j-1)*(sigma4(j)-h2*k34) - theta*(sigma5(j)-h2*k35);
        k43=-1+(mu+alpha+delta1+v(j-1))*(sigma3(j)-h2*k33) - (lamda+v(j-1))*(sigma4(j)-h2*k34);
        k44=(mu+gamma+delta2)*(sigma4(j)-h2*k34) - gamma*(sigma5(j)-h2*k35);
        k45=mu*(sigma5(j)-h2*k35);
        
        sigma1(j-1)=sigma1(j)-(h/6)*(k11+2*k21+2*k31+k41);
        sigma2(j-1)=sigma2(j)-(h/6)*(k12+2*k22+2*k32+k42);
        sigma3(j-1)=sigma3(j)-(h/6)*(k13+2*k23+2*k33+k43);
        sigma4(j-1)=sigma4(j)-(h/6)*(k14+2*k24+2*k34+k44);
        sigma5(j-1)=sigma5(j)-(h/6)*(k15+2*k25+2*k35+k45);
        
    end
    
    count = count +1;
    
    temp1=((sigma2-sigma1)*beta.*S.*E)/(A*N);
    temp2=((sigma3-sigma4).*I)/B;
    temp3=((sigma2-sigma4).*E)/C;
    u11=min(1,max(0,temp1));
    v1=min(1,max(0,temp2));
    w1=min(1,max(0,temp3));
    u=0.5*(u11 + oldu); v=0.5*(v1 + oldv); w=0.5*(w1 + oldw);
    x=[S E I H R]';
    sigma=[sigma1 sigma2 sigma3 sigma4 sigma5]';
    
    test=delta*sum(abs(u))-sum(abs(oldu-u));
    
end

y(1,:)=t;
y(2,:)=S;
y(3,:)=E;
y(4,:)=I;
y(5,:)=H;
y(6,:)=R;
y(7,:)=u;
y(8,:)=sigma1;
y(9,:)=sigma2;
y(10,:)=sigma3;
y(11,:)=sigma4;
y(12,:)=sigma5;


%% plot the evolution of the population E with and without control
fig1 = figure;
hold on
plot(y(1,:),y(3,:),'b',y(1,:),z2(:),'r');
h=legend('$E(t)$ with control','$E(t)$ without control','Location','northeast');
set(h,'Interpreter','latex')
xlabel('Time(days)');  ylabel('E(t)');
title('Exposed with and without control');
%axis([0 15 0 36910558]);
grid on;

%% plot the evolution of the population I with and without control
fig2 = figure;
hold on
plot(y(1,:),y(4,:),'b',y(1,:),z3(:),'r');
h=legend('$I(t)$ with control','$I(t)$ without control','Location','northeast');
set(h,'Interpreter','latex')
xlabel('Time(days)');  ylabel('I(t)');
title('Infected with and without control');
%axis([0 15 0 15000]);
grid on;

%% plot the evolution of the population H with and without control
fig3 = figure;
hold on
plot(y(1,:),y(5,:),'b',y(1,:),z4(:),'r');
h=legend('$H(t)$ with control','$H(t)$ without control','Location','northeast');
set(h,'Interpreter','latex')
xlabel('Time(days)');  ylabel('H(t)');
title('Hospitalized with and without control');
%axis([0 15 0 15000]);
grid on;

%% plot the evolution of the population R with and without control
fig5 = figure;
hold on
plot(y(1,:),y(6,:),'b',y(1,:),z5(:),'r');
h=legend('$R(t)$ with control','$R(t)$ without control','Location','northeast');
set(h,'Interpreter','latex')
xlabel('Time(days)');  ylabel('R(t)');
title('Recovered with and without control');
%axis([0 15 0 500000]);
grid on;

%% plot the evolution of the control u
figure
hold on
plot(y(1,:),u11(:),'color','b');
h=legend('$u(t)$  control','Location','northeast');
set(h,'Interpreter','latex')
xlabel('time (in days)');
ylabel('Control u');
title('(u)');
%axis([0 15 0 1]);
grid on;

%% plot the evolution of the control v
figure
hold on
plot(y(1,:),v1(:),'color','b');
h=legend('$v(t)$  control','Location','northeast');
set(h,'Interpreter','latex')
xlabel('time (in days)');
ylabel('Control v');
title('(v)');
%axis([0 15 0 1]);
grid on;

%% plot the evolution of the control w
figure
hold on
plot(y(1,:),w1(:),'color','b');
h=legend('$w(t)$  control','Location','northeast');
set(h,'Interpreter','latex')
xlabel('time (in days)');
ylabel('Control w');
title('(w)');
%axis([0 15 0 1]);
grid on;

fig3 = figure;

subplot(3,3,1)
plot(y(1,:),z1(:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Susceptible');
hold off
grid
subplot(3,3,2)
plot(y(1,:),z2(:),'c');
hold on
title('(E)');
xlabel('Time(days)');  ylabel('Cumulative number of Exposed');
%axis([0 1500 0 100]);
hold off
grid
subplot(3,3,3)
plot(y(1,:),z3(:),'r');
hold on
title('(I)');
xlabel('Time(days)');  ylabel('Cumulative number of Infected');
hold off
grid
subplot(3,3,4)
plot(y(1,:),z4(:),'y');
hold on
title('(H)');
xlabel('Time(days)');  ylabel('The number of Hospitalized');
hold off
grid
subplot(3,3,5)
plot(y(1,:),z5(:),'g');
hold on
title('(R)');
xlabel('Time(days)');  ylabel('The number of Recovred');
hold off
grid
sgtitle('Population size without control');%comment this line if using previous versions of Matlab 2018b

fig4 = figure;

subplot(3,3,1)
plot(y(1,:),y(2,:),'b');
hold on
title('(S)');
xlabel('Time(days)');  ylabel('Cumulative number of Susceptible');
hold off
grid
subplot(3,3,2)
plot(y(1,:),y(3,:),'c');
hold on
title('(E)');
xlabel('Time(days)');  ylabel('Cumulative number of Exposed');
hold off
grid
subplot(3,3,3)
plot(y(1,:),y(4,:),'r');
hold on
title('(I)');
xlabel('Time(days)');  ylabel('Cumulative number of Infected');
hold off
grid
subplot(3,3,4)
plot(y(1,:),y(5,:),'y');
hold on
title('(H)');
xlabel('Time(days)');  ylabel('The number of Hospitalized');
hold off
grid
subplot(3,3,5)
plot(y(1,:),y(6,:),'g');
hold on
title('(R)');
xlabel('Time(days)');  ylabel('The number of Recovred');
hold off
grid
sgtitle('Population size with control');%comment this line if using previous versions of Matlab 2018b

