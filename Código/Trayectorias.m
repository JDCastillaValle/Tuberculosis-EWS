% Se pueden recorrer los estados estables de salud o los de enefermedad
% con un intervalo de valores equidistantes para el parámetro de fagocitosis. 
% Se van oteniendo los estados estables representantes de salud o enfermedad, 
% según se haya elegido. Estos estados estables se van graficando en la figura 
% generada de "Figueas_de_estados.m" y se utilizan como punto inicial de las trayectorias
% estocásticas, generadas posteriormente en este mismo código, usando el algoritmo de
% Euler-Maruyama. Se obtiene de las trayectorias las medidas estadísticas de 
% varianza, autocorrelación, asimetría, flikcering y varianza de la trayectoria
% hasta antes del flickering.


%%
tic
clear 
Lyapunov

estado= 2; %1 repressenta estados sanos, 2 para enfermedad
if estado == 2
    disp('Recorriendo en enferemdad')
else
    disp('Recorriendo en salud')
end
epsilon = 10.^linspace(-8.5,-6.52579,10);
bif_derecha = 2.9850e-07;
bif_izquierda = .2e-8;

%%
% Fijos
k1=0.00904808526210026; % muerte natural de M
k3=0.0132808065832032; % muerte natural de Mf
k8=2.40598894920150; % tasa de crecimiento de T
k12=2.19806387770886; % tasa de crecimiento de Tf
k15=10384271.4493836; % capacidad de carga de T y de Tf
k16=5.46826266888095; %numero promedio de Tf fagocitada por Mf
k17=18.5909776955088; %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis

%Modulables, elegir una fase (actualmente esta en fase 2)

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

%% Arreglos en donde se guardarán las medidas de cada trayectoria

varianza =[];
varianza_sf = []; %Varianza de la trayectoria hasta el primer flickering
skew=[];
flickering=[];
index=[];

%%
point_ini_cond = 0;
for i = 1:length(epsilon)
    k2=  epsilon(i);% fagocitosis cerquita de la bif derecha

%%
%sacamos los estados estables del valor de k2
syms M_t Mf_t T_t Tf_t

dydt1= 0==Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dydt2= 0==M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dydt3= 0==k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dydt4= 0==k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

equations = [dydt1 dydt2 dydt3 dydt4];
vars=[M_t Mf_t T_t Tf_t];

rango = [NaN NaN; NaN NaN;NaN NaN; NaN NaN];
sol = vpasolve(equations, vars, rango);
ic_stable = [];
ic_unstable = [];

dM=Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dMf=M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dT=k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dTf=k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

J=jacobian([dM dMf dT dTf], vars);

for sol_num = 1:length(sol.M_t)
if isreal([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])
    if sum(double([sol.M_t(sol_num), sol.Mf_t(sol_num), sol.T_t(sol_num), sol.Tf_t(sol_num)]>=0))==4
        positive_ok=1
        Jeval=subs(J, vars, [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
        eigenvals=eig(Jeval);
    if sum(double(eigenvals<0)) == 4
        stable_ok = 1;
        %disp([sol.M_t(sol_num), sol.Mf_t(sol_num), sol.T_t(sol_num), sol.Tf_t(sol_num)])
        if estado == 1
            ic_stable_1 = [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)];
        else
            ic_stable_2 = [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)];
        end
        ic_stable = [ic_stable; sol.M_t(sol_num), sol.Mf_t(sol_num), sol.T_t(sol_num), sol.Tf_t(sol_num)];
        %disp(size(ic_stable))
    else
        stable_ok = 0
        ic_unstable = [ic_unstable; sol.M_t(sol_num), sol.Mf_t(sol_num), sol.T_t(sol_num), sol.Tf_t(sol_num)];
    end
    else
        disp("eigenvals < 0")
    end
else
    disp("no real?")
end
end


if k2 == bif_derecha
    color = 'k';
elseif estado== 1
    color = 'yellow';
else
    color = 'magenta';
end
cond = '';
if isempty(ic_stable)
    %o_O
else
    try
        y = ic_stable(estado,:);
    catch
        y = ic_stable(1,:);
    end

    disp(y)
    
    if y(1)~=0 && estado == 2
        cond = 'solo estado saludable';
    elseif y(1) == 0 && estado == 1
        cond = 'solo estado de enfermedad';
    else
        cond = 'existe biestabilidad';
    end
end

if contains(cond, 'existe')
    %los ploteamos en la gráfica de estados estables positivos
    figure(1)
    hold on
    subplot(2,2,1)
    scatter(log10(k2),y(1), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)

    subplot(2,2,2)
    scatter(log10(k2),y(2), 'filled', 'd', ...
            'MarkerEdgeColor',color, 'MarkerFaceColor', color)

    subplot(2,2,3)
    scatter(log10(k2),y(3), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor', color)

    subplot(2,2,4)
    scatter(log10(k2),y(4), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor', color)
    point_ini_cond = point_ini_cond + 1;
if estado == 1
    path='path';
    name='Recorrido_'+string(i);
    set(gcf, 'Renderer', 'painters');
    saveas(gcf, fullfile(path, name), 'svg');
else
    path='path';
    name='Recorrido_'+string(i);
    set(gcf, 'Renderer', 'painters');
    saveas(gcf, fullfile(path, name), 'svg');
end

unstable_for_flicks = [135257097.511546	0.590452446233121 7501325.42287953 1];

%% Lyapunov
%{
syms k_2
   syms M_t Mf_t T_t Tf_t

   vars=[M_t Mf_t T_t Tf_t];

dM=Mf_t*k6-M_t*(T_t/k16)*k_2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dMf=M_t*(1/k16)*T_t*k_2*(1+Mf_t*k13*k17)  -Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dT=k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k_2*(1+Mf_t*k13*k17);
dTf=k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k_2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;


J=jacobian([dM dMf dT dTf], vars);
%%
% Evaluate the jacobian matrix in the solution
Jeval1=subs(J, vars, y);
 
 %%

A = double(subs(Jeval1, k_2, 0)); 
B = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
Q = lyap(A,B)
clear A B J Jeval1
%}
%% Euler-Maruyama
    realizaciones = 1;

    dt = 0.0001; %incremento del tiempo 0.001
    nPeriods = 30000000; %Numero de intervalos en donde simular 1000000

    F = @(t,X)[X(2)*k6-X(1)*(X(3)/k16)*k2*(1+X(2)*k13*k17)-X(1)*k1-X(1)*X(3)*k9; ...dM
               X(1)*(1/k16)*X(3)*k2*(1+X(2)*k13*k17)-X(2)*k3-X(2)*(X(3)*k10+k14)-X(2)*k5-X(2)*k13; ...dMf
               k8*X(3)*(1-X(3)/k15)+X(2)*k13*k16-X(1)*(X(3)/k16)*k2*(1+X(2)*k13*k17); ...dT
               k12*X(4)*(1-X(4)/(1+k15*X(2)))+X(1)*(1/k16)*X(3)*k2*(1+X(2)*k13*k17)-X(2)*X(4)*k4];...dTf

    if estado == 1
        G =@(t,X) [0.01 0    0   0;
                    0  0.01  0   0;
                    0   0   0.01 0;
                    0   0   0   0.01].*X;% Difussion term : ruido aditivo%

    else 
        G=@(t,X)Q.*X; %usando lyapunov, para el caso de enfermedad
    end

    ic_stable = cast(ic_stable, 'double');


    ini_cond = [y(1);
                y(2);
                y(3);
                y(4)];

    obj = sde(F, G,'StartState', ini_cond) ; % dX = F(t,X)dt + G (t,X)dWa 

    %Simulamos
    [S,T] = simulate(obj , nPeriods , 'DeltaTime', dt , 'nTrials', ...
        realizaciones );

    if any(isinf(S))
        varianza(i,1) = NaN;
        varianza(i,2) = NaN;
        varianza(i,3) = NaN;
        varianza(i,4) = NaN;
    
        [autocor{i,1},~] = NaN;
        [autocor{i,2},~] = NaN;
        [autocor{i,3},~] = NaN;
        [autocor{i,4},~] = NaN;
      
        skew(i,1) = NaN;
        skew(i,2) = NaN;
        skew(i,3) = NaN;
        skew(i,4) = NaN;
    

        flickering(i,1) = NaN;
        flickering(i,2) = NaN;
        flickering(i,3) = NaN;
        flickering(i,4) = NaN;

        index{i,1} = NaN;
        index{i,2} = NaN;
        index{i,3} = NaN;
        index{i,4} = NaN;        

        varianza_sf(i,1) = NaN;
        varianza_sf(i,2) = NaN;
        varianza_sf(i,3) = NaN;
        varianza_sf(i,4) = NaN;
    else 
        varianza(i,1) = nanvar(S(:,1));
        varianza(i,2) = nanvar(S(:,2));
        varianza(i,3) = nanvar(S(:,3));
        varianza(i,4) = nanvar(S(:,4));

        [autocor{i,1},lags] = autocorr(S(:,1));
        [autocor{i,2},lags] = autocorr(S(:,2));

        skew(i,1) = skewness(S(:,1));
        skew(i,2) = skewness(S(:,2));
        skew(i,3) = skewness(S(:,3));
        skew(i,4) = skewness(S(:,4));
    
        value_M  = unstable_for_flicks(1);
        value_Mf = unstable_for_flicks(2);
        value_T  = unstable_for_flicks(3);
        value_Tf = unstable_for_flicks(4);

        [flickering(i,1), index{i,1}] = Flickering(S(:,1),value_M);
        [flickering(i,2), index{i,2}] = Flickering(S(:,2),value_Mf);
        [flickering(i,3), index{i,3}] = Flickering(S(:,3), value_T);
        [flickering(i,4), index{i,4}] = Flickering(S(:,4), value_Tf);
    
        M_sf = min(index{i,1});
        Mf_sf= min(index{i,2});
        T_sf = min(index{i,3});
        Tf_sf= min(index{i,4});
    
    if ~isempty("M_sf")
        varianza_sf(i,1) = nanvar(S(1:M_sf ,1));
    else
        varianza_sf(i,1) = varianza(i,1);
    end
    if ~isempty("Mf_sf")
        varianza_sf(i,2) = nanvar(S(1:Mf_sf,2));
    else
        varianza_sf(i,2) = varianza(i,2);
    end
    if ~isempty("T_sf")
        varianza_sf(i,3) = nanvar(S(1:T_sf ,3));
    else
        varianza_sf(i,3) = varianza(i,3);
    end
    if ~isempty("Tf_sf")
        varianza_sf(i,4) = nanvar(S(1:Tf_sf,4));
    else
        varianza_sf(i,4) = nanvarianza(i,4);
    end
    end


%Graficamos las trayectorias
    if estado==1
    
        figure
        sgtitle(strcat('p2 = ', num2str((k2))))%, ' dt = ', num2str(dt)))
        plot (T, S(:,1),'Color','k');%'#EDB120'
        hold on
        line([0,nPeriods*dt], [unstable_for_flicks(1), ...
            unstable_for_flicks(1)], 'Color', 'b', 'LineStyle', '-')
        line([0,nPeriods*dt], [y(1),y(1)], 'Color', 'g', 'LineStyle', '-')
        ylabel('M')
        xlabel('tiempo')
        xlim([0 3000])
        ylim([0, 2.5e9])
        
        path='path';
        name='M'+string(i);
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

        figure
        sgtitle(strcat('p2 = ', num2str((k2))))%, ' dt = ', num2str(dt)))
        plot (T, S(:,2),'Color','k');%'#EDB120'
        hold on
        line([0,nPeriods*dt], [unstable_for_flicks(2), ...
            unstable_for_flicks(2)], 'Color', 'b', 'LineStyle', '-')
        line([0,nPeriods*dt], [y(2),y(2)], 'Color', 'g', 'LineStyle', '-')
        ylabel('Mf')
        xlabel('tiempo')
        xlim([0 3000])
        ylim([0, 2.5])

        path='path';
        name='Mf'+string(i);
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');


        figure
        sgtitle(strcat('p2 = ', num2str((k2))))
        plot (T, S(:,3),'Color','k');%'#EDB120'
        hold on
        line([0,nPeriods*dt], [unstable_for_flicks(3), ...
            unstable_for_flicks(3)], 'Color', 'b', 'LineStyle', '-')
        line([0,nPeriods*dt], [y(3),y(3)], 'Color', 'g', 'LineStyle', '-')
        ylabel('T')
        xlabel('tiempo')
        xlim([0 3000])
        ylim([0,20000])

        path='path';
        name='T'+string(i);
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');


else
    
    figure
    sgtitle(strcat('p2 = ', num2str((k2))))%, ' dt = ', num2str(dt)))
    subplot(2,1,1)
    axis square
    plot (T, S(:,3),'Color','m');
    hold on
    %for silla = [silla1 silla2 silla3]
    %    if sum(min(ic_stable)<silla)==4 && sum(silla<max(ic_stable))==4
    %        line([0,nPeriods*dt], [silla(3),silla(3)], 'Color', 'k', 'LineStyle', '-')
    %    end
    %end
    %line([0,nPeriods*dt], [silla1(3), silla1(3)], 'Color', '#ff8c00', 'LineStyle', '-')
    %line([0,nPeriods*dt], [silla2(3), silla2(3)], 'Color', '#db7093', 'LineStyle', '-')
    line([0,nPeriods*dt], [silla3(3), silla3(3)], 'Color', 'k', 'LineStyle', '-')%#4DBEEE'%#80B3FF
    try
        min(ic_stable(1,3), ic_stable(2,3))<=ic_unstable(i,3) && max(ic_stable(1,3), ic_stable(2,3))>=ic_unstable(i,3);
        line([0,nPeriods*dt], [ic_unstable(i,3),ic_unstable(i,3)], 'Color', 'b', 'LineStyle', '-')
    catch
        disp('no unstable')
    end
    line([0,nPeriods*dt], [y(3),y(3)], 'Color', 'r', 'LineStyle', '--')
    %line([0,1e-3], [ic(2,1),ic(2,1)], 'Color', 'y', 'LineStyle', '-')
    ylabel('T')
    
    xlabel('tiempo')
    xlim([0 3000])
    ylim([0, 15.5e6])
%{
    subplot(2,2,2)
    histogram(S(:,3))
    xlabel('Densidad de probabilidad')
    ylim([0 10e4])
    xlim([0 1.5e7])
%}


    subplot(2,1,2)
    axis square
    plot (T, S(:,4),'Color','m');
    hold on
    %line([0,nPeriods*dt], [silla1(4), silla1(4)], 'Color', '#ff8c00', 'LineStyle', '-')
    line([0,nPeriods*dt], [silla2(4), silla2(4)], 'Color', 'k', ...
        'LineStyle', '-')%#4DBEEE'%#db7093
    %line([0,nPeriods*dt], [silla3(4), silla3(4)], 'Color', '#80B3FF', 'LineStyle', '-')
    %for silla = [silla1 silla2 silla3]
    %    if sum(min(ic_stable)<silla)==4 && sum(silla<max(ic_stable))==4
    %        line([0,nPeriods*dt], [silla(4),silla(4)], 'Color', 'k', 'LineStyle', '-')
    %    end
    %end
    try min(ic_stable(1,4), ...
            ic_stable(2,4))<=ic_unstable(i,4) && max(ic_stable(1,4), ...
            ic_stable(2,4))>=ic_unstable(i,4);
        line([0,nPeriods*dt], [ic_unstable(i,4),ic_unstable(i,4)], ...
            'Color', 'b', 'LineStyle', '-')
    catch
        disp('unstables out of bounds')
    end
    line([0,nPeriods*dt], [y(4),y(4)], 'Color', 'r', 'LineStyle', '--')
    %line([0,1e-3], [ic(2,2),ic(2,2)], 'Color', 'y', 'LineStyle', '-')
    ylabel('Tf')
    xlabel('tiempo')
    xlim([0 3000])
    ylim([0, 10])

    path='C:\Users\José Daniel\Documents\Tesis\Imágenes\Auto_guardadas\Enfermedad\TyTf';
        name='TyTf'+string(i);
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');
%    xlim([ic_stable(l,4)-5,ic_stable(estado,4)+5])
%{
subplot(2,2,4)
    histogram(S(:,4))
    xlabel('Densidad de probabilidad')
    ylim([0 2e5])
    xlim([0 6])
%}  

    %Totales
%{
    Ttot = S(:,3) + S(:,4);
    
    figure
    sgtitle(strcat('k2 = ', num2str(k2), ' dt = ', num2str(dt)))
    plot(T, Ttot)
    line([0,nPeriods*dt], [ic_stable(estado,3) + ic_stable(estado,4), ...
                           ic_stable(estado,3) + ic_stable(estado,4)], ...
                           'Color', 'r', 'LineStyle', '-')
    line([0,nPeriods*dt], [ic_unstable(1,3) + ic_unstable(1,4), ...
                           ic_unstable(1,3) + ic_unstable(1,4)], ...
                           'Color', '#4DBEEE', 'LineStyle', '--')
    line([0,nPeriods*dt], [ic_unstable(2,3) + ic_unstable(2,4), ...
                           ic_unstable(2,3) + ic_unstable(2,4)], ...
                           'Color', '#4DBEEE', 'LineStyle', '--')
    line([0,nPeriods*dt], [ic_unstable(3,3) + ic_unstable(3,4), ...
                           ic_unstable(3,3) + ic_unstable(3,4)], ...
                           'Color', '#4DBEEE', 'LineStyle', '--')
    line([0,nPeriods*dt], [ic_unstable(4,3) + ic_unstable(4,4), ...
                           ic_unstable(4,3) + ic_unstable(4,4)], ...
                           'Color', '#4DBEEE', 'LineStyle', '--')
    ylabel('Ttot')
    xlabel('tiempo')
%}
   end
else
        varianza(i,1) = NaN;
        varianza(i,2) = NaN;
        varianza(i,3) = NaN;
        varianza(i,4) = NaN;
        
        skew(i,1) = NaN;
        skew(i,2) = NaN;
        skew(i,3) = NaN;
        skew(i,4) = NaN;
    

        flickering(i,1) = NaN;
        flickering(i,2) = NaN;
        flickering(i,3) = NaN;
        flickering(i,4) = NaN;

        index{i,1} = NaN;
        index{i,2} = NaN;
        index{i,3} = NaN;
        index{i,4} = NaN;        

        varianza_sf(i,1) = NaN;
        varianza_sf(i,2) = NaN;
        varianza_sf(i,3) = NaN;
        varianza_sf(i,4) = NaN;    
end

try
disp(size(S))
catch 
    disp('saltamos S')
end

clear T
clear cond
end

[row, col] = find(isnan(varianza));
varianza = rmmissing(varianza,'MinNumMissing',4);
skew = rmmissing(skew,'MinNumMissing',4);
flickering = rmmissing(flickering,'MinNumMissing',4);
index = rmmissing(index,'MinNumMissing',4);
if all(row(1:4)==row(1))
    varianza_sf(row(1), :) = [];
end

%% Graficámos las medidas
if estado == 1
figure
hold on
sgtitle('Varianza')
subplot(2,2,1)
plot(log10(epsilon(1:point_ini_cond)), varianza(:,1), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(:,1)), ...
    max(varianza(:,1))], 'Color', 'k', 'LineStyle', '--') ;%pto. de bif
    ylabel('var(M)')
    xlabel('log10(p2)')

subplot(2,2,2)
plot(log10(epsilon(1:point_ini_cond)), varianza(:,2), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(:,2)), ...
    max(varianza(:,2))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(Mf)')
    xlabel('log10(p2)')

subplot(2,2,[3,4])
plot(log10(epsilon(1:point_ini_cond)), varianza(:,3), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(:,3)), ...
    max(varianza(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(T)')
    xlabel('log10(p2)')

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Salud'];
        name='Varianza';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');
%{
subplot(2,2,3)
plot(log10(epsilon(1:point_ini_cond)), varianza(:,3), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(:,3)), ...
    max(varianza(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(T)')
    xlabel('(p2)')

subplot(2,2,4)
plot(log10(epsilon(1:point_ini_cond)), varianza(:,4), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(:,4)), ...
    max(varianza(:,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(Tf)')
    xlabel('(p2)')
%}

    figure
    hold on
sgtitle('Varianza sin flickering')
subplot(2,2,1)
%plot(log10(epsilon(1:point_ini_cond)), varianza_sf(:,1), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza_sf(:,1)), ...
    max(varianza_sf(:,1))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(M)')
    xlabel('log10(p2)')

subplot(2,2,2)
%plot(log10(epsilon(1:point_ini_cond)), varianza_sf(:,2), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza_sf(:,2)), ...
    max(varianza_sf(:,2))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(Mf)')
    xlabel('log10(p2)')

subplot(2,2,3)
%plot(log10(epsilon(1:point_ini_cond)), varianza_sf(:,3), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza_sf(:,3)), ...
    max(varianza_sf(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(T)')
    xlabel('log10(p2)')

subplot(2,2,4)
%plot(log10(epsilon(1:point_ini_cond)), varianza_sf(:,4), '*-r')
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza_sf(:,4)), ...
    max(varianza_sf(:,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(Tf)')
    xlabel('log10(p2)')

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Salud'];
        name='Varianza sin flickering';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

   
figure
hold on
sgtitle('Flickering')
        subplot(2,2,1)
        plot(log10(epsilon(1:point_ini_cond)), flickering(:,1),'*-g')
        line([log10(bif_izquierda), log10(bif_izquierda)], [min(flickering(:,1)), ...
    max(flickering(:,1))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        ylabel('Número de flicks en M)')
        xlabel('log10(p2)')
    
        subplot(2,2,2)
        plot(log10(epsilon(1:point_ini_cond)), flickering(:,2),'*-g')
        line([log10(bif_izquierda), log10(bif_izquierda)], [min(flickering(:,2)), ...
    max(flickering(:,2))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        ylabel('Número de flicks en  Mf')
        xlabel('log10(p2)')
        
        subplot(2,2,3)
        plot(log10(epsilon(1:point_ini_cond)), flickering(:,3),'*-g')
        line([log10(bif_izquierda), log10(bif_izquierda)], [min(flickering(:,3)), ...
    max(flickering(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        ylabel('Número de flicks en  T')
        xlabel('log10(p2)')
        
        subplot(2,2,4)
        plot(log10(epsilon(1:point_ini_cond)), flickering(:,4),'*-g')
        line([log10(bif_izquierda), log10(bif_izquierda)], [min(flickering(:,4)), ...
    max(flickering(:,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        ylabel('Número de flicks en Tf')
        xlabel('log10(p2)')

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Salud'];
        name='Flickering';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

figure
hold on
sgtitle('Asimetría')
    subplot(2,2,1)
    plot(log10(epsilon(1:point_ini_cond)),skew(:,1),'*-b')
    line([log10(bif_izquierda), log10(bif_izquierda)], [min(skew(:,1)), ...
    max(skew(:,1))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('Asimetría de M')
    xlabel('log10(p2)')

    subplot(2,2,2)
    plot(log10(epsilon(1:point_ini_cond)),skew(:,2),'*-b')
    line([log10(bif_izquierda), log10(bif_izquierda)], [min(skew(:,2)), ...
    max(skew(:,2))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('Asimetría de Mf')
    xlabel('log10(p2)')

    subplot(2,2,[3,4])
    plot(log10(epsilon(1:point_ini_cond)),skew(:,3),'*-b')
    line([log10(bif_izquierda), log10(bif_izquierda)], [min(skew(:,3)), ...
    max(skew(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('Asimetría de T')
    xlabel('log10(p2)')  

        path='path';
        name='Asimetría';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');
%{
    subplot(2,2,4)
    plot(log10(epsilon(1:point_ini_cond)),skew(:,4),'*-b')
    line([log10(bif_izquierda), log10(bif_izquierda)], [min(skew(:,4)), ...
    max(skew(:,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('skewness of Tf')
    xlabel('(p2)')
%}
else

    figure
    hold on
sgtitle('Varianza')
subplot(2,1,1)
plot(log10(epsilon(1:point_ini_cond)), varianza(:,3), '*-r')
line([log10(2.98e-7), log10(2.98e-7)], [min(varianza(:,3)), ...
    max(varianza(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(T)')
        xlabel('log10(p2)')

subplot(2,1,2)
plot(log10(epsilon(1:point_ini_cond)), varianza(:,4), '*-r')
line([log10(2.98e-7), log10(2.98e-7)], [min(varianza(:,4)), ...
    max(varianza(:,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(Tf)')
    xlabel('log10(p2)')

        path'path';
        name='Varianza';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

     figure
     hold on
sgtitle('Varianza sin flickering')
subplot(2,1,1)
plot(log10(epsilon(1:point_ini_cond)), varianza_sf(:,3), '*-r')
line([log10(2.98e-7), log10(2.98e-7)], [min(varianza_sf(:,3)), ...
    max(varianza_sf(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(T)')
        xlabel('log10(p2)')

subplot(2,1,2)
plot(log10(epsilon(1:point_ini_cond)), varianza_sf(:,4), '*-r')
line([log10(2.98e-7), log10(2.98e-7)], [min(varianza_sf(:,4)), ...
    max(varianza_sf(:,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('var(Tf)')
    xlabel('log10(p2)')

        path='path';
        name='Varianza sin flickering';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

    figure
    hold on
sgtitle('Flickering')

        subplot(2,1,1)
        plot(log10(epsilon(1:point_ini_cond)), flickering(:,3),'*-g')
        line([log10(2.98e-7), log10(2.98e-7)], [min(flickering(:,3)), ...
    max(flickering(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        ylabel('Número de flicks en T')
            xlabel('log10(p2)')

        subplot(2,1,2)
        plot(log10(epsilon(1:point_ini_cond)), flickering(:,4),'*-g')
        line([log10(2.98e-7), log10(2.98e-7)], [min(flickering(:,4)), ...
    max(flickering(:,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        ylabel('Número de flicks en Tf')
        xlabel('log10(p2)')

        path='path';
        name='Flickering';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');



figure
hold on
sgtitle('Asimetría')
    subplot(2,1,1)
    plot(log10(epsilon(1:point_ini_cond)),skew(:,3),'*-b')
    line([log10(2.98e-7), log10(2.98e-7)], [min(skew(:,3)), ...
    max(skew(:,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('Asimetría de T')
    xlabel('log10(p2)')

    subplot(2,1,2)
    plot(log10(epsilon(1:point_ini_cond)),skew(:,4),'*-b')
    line([log10(2.98e-7), log10(2.98e-7)], [min(skew(:,4)), ...
    max(skew(:,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('Asimetría de Tf')
    xlabel('log10(p2)')
end
        path='path';
        name='Asimetría';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');
toc
