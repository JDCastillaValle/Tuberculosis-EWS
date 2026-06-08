%% Versión completa, tardada

%
close all
clear all
clc

%% Elegir una fase 
			
tic

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

%%
p1=k6;  %  Recruitment of macrophages to the site of infection. PHASE DEPENDENT
p3=k13*k17; % re-release of bacteria for phagocytosis by necrotizing Mf PHASE DEPENDENT (K13 is phase dependent)
p4=k1;   % natural death of M ---- CONSTANT
p5=k9;   % death of M by T PHASE DEPENDENT
p6=k3+k14+k5+k13; % linear death (various forms) of Mf - PHASE DEPENDENT (K13 AND K14)
p7=k10; %  death of Mf by T PHASE DEPENDENT
p8=k8; % growth rate of T constant
p9=k15; % carriyng capacity of T and of Tf% CONSTANT
p10=k13*k16; % release of T through necrosis PHASE DEPENDENT
p11=k12; % growth rate of Tf constant
p12=k4;

%% 


k2=[0, 1e-10, .75e-10, .5e-10, .25e-10, ...
        1e-9,  .75e-9, .5e-9, .25e-9,  ...
        1e-8,  .75e-8, .5e-8, .25e-8,  ...
        1e-7, .75e-7, .5e-7, .25e-7,  ...
        1e-6, .75e-6, .5e-6, .25e-6,  ...
        1e-5, .75e-5, .5e-5, .25e-5,  ...
        1e-4, .75e-4, .5e-4, .25e-4];

for i=1:length(k2)

   p2=k2(i)/k16; % phagocytosis rate PHASE DEPENDENT

    syms M_t Mf_t T_t Tf_t


      % recruit       % phagocytosis     % natural death    % death by T
dydt1= 0==Mf_t*p1 - M_t*T_t*p2*(1+Mf_t*p3) -M_t*p4           - M_t*T_t*p5;
                  % phagocytosis       % linear death    % death by T
dydt2= 0==          M_t*T_t*p2*(1+Mf_t*p3) -Mf_t*p6          -Mf_t*T_t*p7;

    % logistic growth   % release through necrosis     % phagocytosis   
dydt3= 0== p8*T_t*(1-T_t/p9)       + Mf_t*p10                  -M_t*T_t*p2*(1+Mf_t*p3);
      % logistic growth             % phagocytosis         % death of Tf by Mf
dydt4= 0==p11*Tf_t*(1-Tf_t/(1+p9*Mf_t))+ M_t*T_t*p2*(1+Mf_t*p3) -Mf_t*Tf_t*p12;


    equations = [dydt1 dydt2 dydt3 dydt4];
    vars=[M_t Mf_t T_t Tf_t];

    range = [NaN NaN; NaN NaN;NaN NaN; NaN NaN];
    sol = vpasolve(equations, vars, range);


       % recruit       % phagocytosis     % natural death    % death by T
dM = Mf_t*p1 - M_t*T_t*p2*(1+Mf_t*p3) -M_t*p4           - M_t*T_t*p5;
                       % phagocytosis     % linear death     % death by T
dMf=           M_t*T_t*p2*(1+Mf_t*p3) -Mf_t*p6          -Mf_t*T_t*p7;

    % logistic growth   % release through necrosis     % phagocytosis   
dT= p8*T_t*(1-T_t/p9)       + Mf_t*p10                  -M_t*T_t*p2*(1+Mf_t*p3);
      % logistic growth             % phagocytosis         % death of Tf by Mf
dTf=p11*Tf_t*(1-Tf_t/(1+p9*Mf_t))+ M_t*T_t*p2*(1+Mf_t*p3) -Mf_t*Tf_t*p12;


    J=jacobian([dM dMf dT dTf], vars);
%% Verificar qué estados se desean graficar
%estables positivos reales
    stable_positive_real_solution_matrix=[];

% Des-comentar para graficar diferentes estados al real, positivo y estable
% tambíen descomentar más abajo la sección de los estados deseados

%inestables positivos reales
    unstable_positive_real_solution_matrix=[];

%estables negativos reales
    stable_non_positive_real_solution_matrix=[];

%inestables positivos reales
    unstable_non_positive_real_solution_matrix=[];

    conter_stable_pos_real=0;
    % Se discrimina el tipo de estado
    for sol_num=1:1:length(sol.M_t)

        if isreal([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])
            real_ok=1;


            if sum(double(([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])>=0))==4
                positive_ok=1;
            else
                positive_ok=0;
            end


            Jeval=subs(J, vars, [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
            eigenvals=eig(Jeval) ;
            if (sum(double(eigenvals)<0))==4

                stable_ok=1;
                sol_num;

            else
                stable_ok=0;
            end


            
            % Se guardan los diferentes estados estables positivos reales
            if real_ok==1 && positive_ok==1 && stable_ok==1
                stable_positive_real_solution_matrix=[stable_positive_real_solution_matrix; 
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
                
                M_sol_stable(i)={stable_positive_real_solution_matrix(:,1)};
                Mf_sol_stable(i)={stable_positive_real_solution_matrix(:,2)};
                T_sol_stable(i)={stable_positive_real_solution_matrix(:,3)};
                Tf_sol_stable(i)={stable_positive_real_solution_matrix(:,4)};
                
                
                conter_stable_pos_real=conter_stable_pos_real+1;
%Des-comentar para graficar estados inestables positivos y reales

            elseif real_ok==1 && positive_ok==1 && stable_ok==0
                unstable_positive_real_solution_matrix=[unstable_positive_real_solution_matrix;
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
                M_sol_unstable(i)={unstable_positive_real_solution_matrix(:,1)};
                Mf_sol_unstable(i)={unstable_positive_real_solution_matrix(:,2)};
                T_sol_unstable(i)={unstable_positive_real_solution_matrix(:,3)};
                Tf_sol_unstable(i)={unstable_positive_real_solution_matrix(:,4)};
                
                simbolo='d';
                color_afuera='#4DBEEE';
                color_adentro= '#4DBEEE';
                size=25; 
                figNum=2;
            figure(figNum)
            sgtitle('Non-negative real unstable states')
            subplot(2,2,1)
            scatter(log10(p2), (double(sol.M_t(sol_num))), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('(M_t)')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-10, 3e9], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-10, 3e9], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,2)
            scatter(log10(p2),(double(sol.Mf_t(sol_num))), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('(Mf_t)')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-1, 30], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-1, 30], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,3)
            scatter(log10(p2),double(sol.T_t(sol_num)./k15), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('T_t/k15')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-0.1, 1], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-0.1, 1], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,4)
            scatter(log10(p2),double(sol.Tf_t(sol_num)), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('Tf_t')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-1, 40], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-1, 40], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis 

%% Des-comentar para graficar estados estables, negativos y reales

            elseif   real_ok==1 && positive_ok==0 && stable_ok==1 % negative elements, which can be stable or not stable
                simbolo='s';
                color_afuera='b';
                color_adentro= 'b';
                size=50;
                stable_non_positive_real_solution_matrix=[stable_non_positive_real_solution_matrix;
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
                figNum=3;
                
                figure(figNum)
            sgtitle('Negative real stable states')
            subplot(2,2,1)
            scatter(log10(p2), (double(sol.M_t(sol_num))), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('(M_t)')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-10, 3e9], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-10, 3e9], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,2)
            scatter(log10(p2),(double(sol.Mf_t(sol_num))), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('(Mf_t)')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-1, 30], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-1, 30], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,3)
            scatter(log10(p2),double(sol.T_t(sol_num)./k15), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('T_t/k15')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-0.1, 1], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-0.1, 1], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,4)
            scatter(log10(p2),double(sol.Tf_t(sol_num)), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('Tf_t')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-1, 40], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-1, 40], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosi

%% Des-comentar para graficar estados inestables, negativos y reales

            elseif   real_ok==1 && positive_ok==0 && stable_ok==0 % negative elements, which can be stable or not stable
                simbolo='*';
                color_afuera='m';
                color_adentro= 'm';
                size=25;
                unstable_non_positive_real_solution_matrix=[unstable_non_positive_real_solution_matrix;
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
                figNum=4;

                figure(figNum)
            sgtitle('Negative real unstable states')
            subplot(2,2,1)
            scatter(log10(p2), (double(sol.M_t(sol_num))), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('(M_t)')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-10, 3e9], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-10, 3e9], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,2)
            scatter(log10(p2),(double(sol.Mf_t(sol_num))), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('(Mf_t)')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-1, 30], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-1, 30], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,3)
            scatter(log10(p2),double(sol.T_t(sol_num)./k15), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('T_t/k15')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-0.1, 1], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-0.1, 1], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis

            subplot(2,2,4)
            scatter(log10(p2),double(sol.Tf_t(sol_num)), size,'MarkerEdgeColor',color_afuera,...
                'MarkerFaceColor',color_adentro,...
                'LineWidth',1.5, 'Marker', simbolo )
            hold on
            ylabel('Tf_t')
            xlabel('log10(k2)')
            line([log10(2.98e-7), log10(2.98e-7)], [-1, 40], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            line([log10(.2e-8), log10(.2e-8)], [-1, 40], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
            
%}
            
            end % end of clasification of states

        else
            real_ok=0;
        end
        % filtering out the complex solutions


    end % end looping over all the solutions
 end % end looping over k2 values
 
%% Categorizamos estados estables positivos reales

for i = 1:length(M_sol_stable)
    X = sort(M_sol_stable{:,i});
    if length(X) == 1
        M_sol_stable_2(i) = X(1);
        M_sol_stable_1(i) = NaN;
    elseif length(X) == 2
        M_sol_stable_2(i) = X(1);
        M_sol_stable_1(i) = X(2);
    else
        warning('Hay mas de 2 soluciones en M')
    end
end

for i = 1:length(Mf_sol_stable)
    X= sort(Mf_sol_stable{:,i});
    if length(X) == 1
        Mf_sol_stable_2(i) = X(1);
        Mf_sol_stable_1(i) = NaN;
    elseif length(X) == 2
        Mf_sol_stable_2(i) = X(1);
        Mf_sol_stable_1(i) = X(2);
    else
        warning('Hay mas de 2 soluciones en Mf')
    end
end

T_sol_stable_2(1) = cell2mat(T_sol_stable(1));
T_sol_stable_1(1) = NaN;
for i = 2:length(T_sol_stable)
    X = sort(T_sol_stable{:,i});
    if length(X) == 1
        if abs(T_sol_stable_1(i-1) - X(1)) < 1000
            T_sol_stable_1(i) = X(1);
            T_sol_stable_2(i) = NaN;
        else
            T_sol_stable_2(i) = X(1);
            T_sol_stable_1(i) = NaN;
        end

    elseif length(X) == 2
        T_sol_stable_1(i) = X(1);
        T_sol_stable_2(i) = X(2);
    else
        warning('Hay mas de 2 soluciones en T')
    end
end

Tf_sol_stable_2(1) = cell2mat(Tf_sol_stable(1));
Tf_sol_stable_1(1) = NaN;
for i = 2:length(Tf_sol_stable)
    X = sort(Tf_sol_stable{:,i});
    if length(X) == 1
        if abs(Tf_sol_stable_1(i-1) - X(1)) < 1000
            Tf_sol_stable_1(i) = X(1);
            Tf_sol_stable_2(i) = NaN;
        else
            Tf_sol_stable_2(i) = X(1);
            Tf_sol_stable_1(i) = NaN;
        end

    elseif length(X) == 2
        Tf_sol_stable_1(i) = X(1);
        Tf_sol_stable_2(i) = X(2);
    else
        warning('Hay mas de 2 soluciones en Tf')
    end
end
%% Graficamos los estados estables positivos reales
figure(1)
sgtitle('Estados estables no-negativos')
subplot(2,2,1)
hold on
scatter(log10(k2), M_sol_stable_1, 'filled', 'o', 'g')
scatter(log10(k2), M_sol_stable_2, 'filled', 'o', 'r')
line([log10(2.98e-7), log10(2.98e-7)], [0, 2.5e9], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
line([log10(.2e-8), log10(.2e-8)], [0, 2.5e9], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
text(-9,-29e7,'p_2^-')
text(-6.8,-29e7,'p_2^+')
ylabel('M_t')
xlabel('log10(p2)')

subplot(2,2,2)
hold on
scatter(log10(k2), Mf_sol_stable_1, 'filled', 'o', 'g')
scatter(log10(k2), Mf_sol_stable_2, 'filled', 'o', 'r')
line([log10(2.98e-7), log10(2.98e-7)], [0, 2.5], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
line([log10(.2e-8), log10(.2e-8)], [0, 2.5], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
text(-9,-0.3,'p_2^-')
text(-6.8,-0.3,'p_2^+')
ylabel('Mf_t')
xlabel('log10(p2)')

subplot(2,2,3)
hold on
scatter(log10(k2), T_sol_stable_1, 'filled', 'o', 'g')
scatter(log10(k2), T_sol_stable_2, 'filled', 'o', 'r')
line([log10(2.98e-7), log10(2.98e-7)], [0, 11e6], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
line([log10(.2e-8), log10(.2e-8)], [0, 11e6], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
%line([log10(k2(9)), log10(k2(9))], [0, 11e6], 'Color', 'b', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
text(-9,-9e5,'p_2^-')
text(-6.8,-9e5,'p_2^+')
ylabel('T_t')
xlabel('log10(p2)')

subplot(2,2,4)
hold on
scatter(log10(k2), Tf_sol_stable_1, 'filled', 'o', 'g')
scatter(log10(k2), Tf_sol_stable_2, 'filled', 'o', 'r')
line([log10(2.98e-7), log10(2.98e-7)], [0, 1], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
line([log10(.2e-8), log10(.2e-8)], [0, 1], 'Color', 'k', 'LineStyle', '--') ;%1.01382645792285e-10*184.71337580; % fagocitosis
text(-9,-0.09,'p_2^-')
text(-6.8,-0.09,'p_2^+')
ylabel('Tf_t')
xlabel('log10(p2)')

%% Guardamos la imágen en tipo svg
        path='C:\Users\José Daniel\Documents\Tesis\Imágenes\Auto_guardadas';
        name='Figura_de_estados';
        % Set renderer to painters for proper vector graphics
        set(gcf, 'Renderer', 'painters');

        % Save the figure with the full path
        saveas(gcf, fullfile(path, name), 'svg');

%% En caso de haber graficado estados no-negativos reales e inestables
if ishandle(2)
    figure(2)
        subplot(2,2,1)
            text(-9,-0.2e9,'p_2^-')
            text(-6.8,-0.2e9,'p_2^+')
        subplot(2,2,2)
            text(-9,-2.5,'p_2^-')
            text(-6.8,-2.5,'p_2^+')
        subplot(2,2,3)
            text(-9,-0.17,'p_2^-')
            text(-6.8,-0.17,'p_2^+')
        subplot(2,2,4)
            text(-9,-3,'p_2^-')
            text(-6.8,-3,'p_2^+')
end   

%% Des-comentar para graficar estados inestables, negativos y reales




toc
