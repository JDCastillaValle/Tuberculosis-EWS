% Use of Lyapunov equation for significant noise

%% fixed parameters


% Fijos
k1=0.00904808526210026; % muerte natural de M
k3=0.0132808065832032; % muerte natural de Mf
k8=2.40598894920150; % tasa de crecimiento de T
k12=2.19806387770886; % tasa de crecimiento de Tf
k15=10384271.4493836; % capacidad de carga de T y de Tf
k16=5.46826266888095; %numero promedio de Tf fagocitada por Mf
k17=18.5909776955088; %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis

%Modulables, elegir una fase (actualmente esta en fase 1)

%Fase 1

%{

k2=1.01382645792285e-10; % fagocitosis
k4=1099.11111656674; %muerte de Tf por Mf
k5=2610.49082476506; % apoptosis
k6=1835831.53903044; % reclutamiento de M por Mf
k9=1.27012642169026e-10; % muerte de M por T
k10=0.0428656832209084; % muerte de Mf por T
k13=4.02705861695495e-05; % necrosis
k14=306.222081946196; % muerte de Mf por Tf
%}
%Fase 2

%k2= 2.9850e-07;%1.01382645792285e-10*184.71337580; % fagocitosis

k4=1099.11111656674*313.71294168; %muerte de Tf por Mf
k5=2610.49082476506*0.35602934; % apoptosis
k6=1835831.53903044*6.60758304; % reclutamiento de M por Mf
k9=1.27012642169026e-10*7.20104815; % muerte de M por T
k10=0.0428656832209084*26.38000000; % muerte de Mf por T
k13=4.02705861695495e-05*1000; % necrosis
k14=306.222081946196*8.74477588; % muerte de Mf por Tf


%Fase 3
%{
k2=1.01382645792285e-10*184.71337580*0.003; % fagocitosis
k4=1099.11111656674*313.71294168*0.0005; %muerte de Tf por Mf
k5=2610.49082476506*0.35602934*0.005; % apoptosis
k6=1835831.53903044*6.60758304*0.0001; % reclutamiento de M por Mf
k9=1.27012642169026e-10*7.20104815*0.147; % muerte de M por T
k10=0.0428656832209084*26.38000000*0.189; % muerte de Mf por T
k13=4.02705861695495e-05*1000*102.593; % necrosis
k14=306.222081946196*8.74477588*0.302; % muerte de Mf por Tf

%}

%% compute the solutions 
 k2=0;
  syms M_t Mf_t T_t Tf_t


dydt1= 0==Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dydt2= 0==M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dydt3= 0==k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dydt4= 0==k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

equations = [dydt1 dydt2 dydt3 dydt4];
vars=[M_t Mf_t T_t Tf_t];

range = [NaN NaN; NaN NaN;NaN NaN; NaN NaN];
sol = vpasolve(equations, vars, range);

%}
solution_matrix=[sol.M_t sol.Mf_t sol.T_t sol.Tf_t]

%% we have 4 positive real solutions
sol_1=[0,0,0,0]; % matrix independent of k2, always unstable

sol_2=[0,0,0,1]; % % matrix independent of k2, always unstable

sol_3=[0,0,k15,0]; %% two of the eigenvalues is dependent of k2 but one of them is always positive the other one always negative 

sol_4=[0,0,k15,1]; %%% two of the eigenvalues is dependent of k2 ; the remaining two are always negative

%% to asses the stability we compute the jacobian matrix

syms k2
   syms M_t Mf_t T_t Tf_t

   vars=[M_t Mf_t T_t Tf_t];

dM=Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dMf=M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dT=k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dTf=k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;


J=jacobian([dM dMf dT dTf], vars);

%% now we will check each of the solutions 
 Jeval1=subs(J, vars, sol_1)
 Jeval1=double(Jeval1)
 eigenvals1=(eig(Jeval1)) 
%%
 Jeval2=subs(J, vars, sol_2)
 Jeval2=double(Jeval2)
 eigenvals2=(eig(Jeval2))
 %%
  Jeval3=subs(J, vars, sol_3)
%  Jeval3=subs(Jeval3, k2, 0)
%  Jeval3=double(Jeval3)
 eigenvals3=(eig(Jeval3))
%%
%eigenvals3=subs(eigenvals3, k2, 0)
%%
% for k2=0  - + - -
double(subs(eigenvals3(3), k2, 0)) %disque
%double(subs(eigenvals3(4), k2, 0))

%%
  Jeval4=subs(J, vars, sol_4)
 eigenvals4=(eig(Jeval4))
 %%
double(subs(eigenvals4(3),k2,0)) %negative
double(subs(eigenvals4(4),k2,0)) %negative
%%
double(solve( eigenvals4(3)==0, k2)) % no solution, always negative
%double(solve( eigenvals4(4)==0, k2)) % change for 2.9850e-07
%% Ecuación de Lyapunov AQ+QA^T+BB^T=0
A = double(subs(Jeval4, k2, 0)); 
B = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
Q = lyap(A,B)
clear A B J Jeval1 Jeval2 Jeval3 