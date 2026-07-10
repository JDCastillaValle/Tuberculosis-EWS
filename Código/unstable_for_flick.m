%% Solo necesitas elegir una fase y correr el codigo, los puntos de equilibrio estables se imprimiran en la consola
			
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

k2=1.01382645792285e-10*184.71337580; % fagocitosis
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

%%
  syms M_t Mf_t T_t Tf_t


dydt1= 0==Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dydt2= 0==M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dydt3= 0==k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dydt4= 0==k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

equations = [dydt1 dydt2 dydt3 dydt4];
vars=[M_t Mf_t T_t Tf_t];

range = [NaN NaN; NaN NaN;NaN NaN; NaN NaN];
sol = vpasolve(equations, vars, range);

%%
dM=Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dMf=M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dT=k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dTf=k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

J=jacobian([dM dMf dT dTf], vars);

stable_positive_real_solution_matrix=[];
unstable_positive_real_solution_matrix=[];

for sol_num=1:1:length(sol.M_t)

%[sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]
if isreal([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])
     if sum(double(([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])>=0))==4
            Jeval=subs(J, vars, [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
        eigenvals=eig(Jeval) ;
           if (sum(double(eigenvals)<0))==4
           disp('solution is stable :)');
           sol_num
            stable_positive_real_solution_matrix=[stable_positive_real_solution_matrix; 
                          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])]
                       
       else
            disp('solution is unstable :(');
            unstable_positive_real_solution_matrix=[unstable_positive_real_solution_matrix; 
               double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
        end
       
           [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)];

     end
end
end
disp(unstable_positive_real_solution_matrix)
unstable_for_flickering = unstable_positive_real_solution_matrix(4,:);
