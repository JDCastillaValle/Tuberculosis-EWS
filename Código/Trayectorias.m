% Se pueden recorrer los estados estables de salud o los de enefermedad
% con un intervalo de valores equidistantes para el parámetro de fagocitosis. 
% Se van oteniendo los estados estables representantes de salud o enfermedad, 
% según se haya elegido. Estos estados estables se van graficando en la figura 
% generada de "Figueas_de_estados.m" y se utilizan como punto inicial de las trayectorias
% estocásticas, generadas posteriormente en este mismo código, usando el algoritmo de
% Euler-Maruyama. Se obtiene de las trayectorias las medidas estadísticas de 
% varianza, autocorrelación, asimetría, flikcering y varianza de la trayectoria
% hasta antes del flickering.



tic
Figuras_de_estados

clear 
Lyapunov

estado= 1; %1 repressenta estados sanos, 2 para enfermedad
if estado == 2
    disp('Recorriendo en enferemdad')
else
    disp('Recorriendo en salud')
end
epsilon = 10.^linspace(-8.5,-6.52579,10);
bif_derecha = 2.9850e-07;
bif_izquierda = 0.2e-8;

%%

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


%% Arreglos en donde se guardarán las medidas de cada trayectoria
varianza = NaN(length(epsilon),4);
varianza_sf = NaN(length(epsilon),4); %Varianza de la trayectoria hasta el primer flickering
skew = NaN(length(epsilon),4);
flickering = NaN(length(epsilon),4);
%index = NaN(length(epsilon),4);
color_flick = repmat([6 64 43], length(epsilon),1);
unstables_for_flick = repmat([135257097.511546	0.590452446233121 7501325.42287953 24.6114869148956],length(epsilon),1);
%%
point_ini_cond = 0;
epsilon_index = [];
for i = 1:length(epsilon)
   p2=  epsilon(i);% fagocitosis cerquita de la bif derecha
   k2 = p2;
   % phagocytosis rate PHASE DEPENDENT


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
negative_stables = [];
negative_unstables = [];
eigenvals=[];

dM=Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dMf=M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dT=k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dTf=k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;


J=jacobian([dM dMf dT dTf], vars);
conter_stable_pos_real = 0;

    for sol_num=1:1:length(sol.M_t)

        %[sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]
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



            % Guardamos las soluciones reales según su estabilidad y signo

            if real_ok==1 && positive_ok==1 && stable_ok==1
                ic_stable = [ic_stable; 
                    sol.M_t(sol_num), sol.Mf_t(sol_num), sol.T_t(sol_num), sol.Tf_t(sol_num)];
                

                conter_stable_pos_real=conter_stable_pos_real+1;
            elseif real_ok==1 && positive_ok==1 && stable_ok==0
                ic_unstable=[ic_unstable;
          double([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
                if sum(double(([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])>=0))==4
                        unstables_for_flick(i,:) = [
                double([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
                        color_flick(i,:) = [28, 255, 24];
                end
            elseif   real_ok==1 && positive_ok==0 && stable_ok==1 % negative elements, which can be stable or not stable
                disp('Negativos estables')
                negative_stables = [negative_stables;
                    double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];

            elseif   real_ok==1 && positive_ok==0 && stable_ok==0 % negative elements, which can be stable or not stable
                disp('Negativos inestables')
                negative_unstables = [negative_unstables;
                    double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];

            end

        else
            real_ok=0;
        end
        % filtering out the complex solutions


    end % end looping over all the solutions


if k2 == bif_derecha
    color = 'k';
else
    color = 'yellow';
end


cond = '';


if isempty(ic_stable)
    %O_o
else
    try
        y = ic_stable(estado,:);
    catch
        y = ic_stable(1,:);
    end
    
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
    epsilon_index = [epsilon_index; i];
    point_ini_cond = point_ini_cond + 1;

    screensize = get(0, "ScreenSize");
 
if estado == 1
    figure(1)
    hold on
    subplot(2,2,1)
    scatter(log10(k2),log10(y(1)), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)
    text(log10(k2), log10(y(1)), string(10-point_ini_cond))

    subplot(2,2,2)
    scatter(log10(k2),log10(y(2)), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)
    text(log10(k2), log10(y(2)), string(10-point_ini_cond))

    subplot(2,2,3)
    scatter(log10(k2),y(3), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)
    text(log10(k2), y(3), string(10-point_ini_cond))

    subplot(2,2,4)
    scatter(log10(k2),y(4), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)
    text(log10(k2), y(4), string(10-point_ini_cond))

    path='C:\Users\José Daniel\Documents\Tesis\Imágenes\Auto_guardadas\Salud';
    name='Recorrido_'+string(i);
    set(gcf, 'Renderer', 'painters', 'Position', screensize);
    saveas(gcf, fullfile(path, name), 'svg');
 
else

    figure(1)
    hold on
    subplot(2,2,1)
    scatter(log10(k2),y(1), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)
    text(log10(k2), y(1), string(point_ini_cond))

    subplot(2,2,2)
    scatter(log10(k2),y(2), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)
    text(log10(k2), y(2), string(point_ini_cond))

    subplot(2,2,3)
    scatter(log10(k2),y(3), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)
    text(log10(k2), y(3), string(point_ini_cond))

    subplot(2,2,4)
    scatter(log10(k2),y(4), 'filled', 'd', ...
            'MarkerEdgeColor', color, 'MarkerFaceColor',color)
    text(log10(k2), y(4), string(point_ini_cond))

    path='C:\Users\José Daniel\Documents\Tesis\Imágenes\Auto_guardadas\Enfermedad';
    name='Recorrido_'+string(i);
    set(gcf, 'Renderer', 'painters', 'Position', screensize);
    saveas(gcf, fullfile(path, name), 'svg');

end

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

A = double(subs(Jeval1, k_2, epsilon(i))); 
B = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
Q = lyap(A,B);

clear A B J Jeval1
%}
%% Euler-Maruyama
    realizaciones = 1;

    dt = 0.0001;%1/max(abs(double(eigenvals))); %incremento del tiempo 0.0001
    nPeriods = 30000000;%floor(max(abs(double(eigenvals)))); %Numero de intervalos en donde simular 1000000

    F = @(t,X)[X(2)*k6-X(1)*(X(3)/k16)*k2*(1+X(2)*k13*k17)-X(1)*k1-X(1)*X(3)*k9; ...dM
               X(1)*(1/k16)*X(3)*k2*(1+X(2)*k13*k17)-X(2)*k3-X(2)*(X(3)*k10+k14)-X(2)*k5-X(2)*k13; ...dMf
               k8*X(3)*(1-X(3)/k15)+X(2)*k13*k16-X(1)*(X(3)/k16)*k2*(1+X(2)*k13*k17); ...dT
               k12*X(4)*(1-X(4)/(1+k15*X(2)))+X(1)*(1/k16)*X(3)*k2*(1+X(2)*k13*k17)-X(2)*X(4)*k4];...dTf

    if estado == 1
        G =@(t,X) [ 0.01 0   0   0;
                    0   0.01 0   0;
                    0    0  0.01 0;
                    0    0   0  0.01].*X;%@(t,X)Q.*X;%@(t,X) [55.7136 0   0     0  ;
            %0  0 0     0  ;
            %0    0  0.2078   0  ;
            %0    0   0   1.0607].*X ; % Difussion term : ruido aditivo%

    else 
        G=@(t,X)Q.*X;
    end

    ic_stable = cast(ic_stable, 'double');


    ini_cond = [y(1);
                y(2);
                y(3);
                y(4)];

    %ini_cond = double(ini_cond);

    obj = sde(F, G,'StartState', ini_cond) ; % dX = F(t,X)dt + G (t,X)dWa 

    %Simulamos
    [S,T] = simulate(obj , nPeriods , 'DeltaTime', dt , 'nTrials', ...
        realizaciones );

 %%
        varianza(i,1) = var(S(isfinite(S(:,1)),1));
        varianza(i,2) = var(S(isfinite(S(:,2)),2));
        varianza(i,3) = var(S(isfinite(S(:,3)),3));
        varianza(i,4) = var(S(isfinite(S(:,4)),4));

        try
    [autocor{i,1},lags] = autocorr(S(isfinite(S(:,1)),1), 'NumLags', 2);
    [autocor{i,2},lags] = autocorr(S(isfinite(S(:,2)),2), 'NumLags', 2);
    [autocor{i,3},lags] = autocorr(S(isfinite(S(:,3)),3), 'NumLags', 2);
    [autocor{i,4},lags] = autocorr(S(isfinite(S(:,4)),4), 'NumLags', 2);
        catch
            numLags = 2;
            %if estado == 1
                for v = 1:1:4
                    x = S(isfinite(S(:,v)),v);
                    N = numel(x);
                    xc = x - mean(x);
                    acf = zeros(numLags+1,1);
                        for k = 0:numLags
                            acf(k+1) = sum(xc(1+k:N) .* xc(1:N-k));
                        end
                        acf = acf / acf(1);
                    autocor{i,v} = acf(2);
                end
                
         end

        skew(i,1) = skewness(S(isfinite(S(:,1)),1));
        skew(i,2) = skewness(S(isfinite(S(:,2)),2));
        skew(i,3) = skewness(S(isfinite(S(:,3)),3));
        skew(i,4) = skewness(S(isfinite(S(:,4)),4));

        value_M  = unstables_for_flick(i,1);
        value_Mf = unstables_for_flick(i,2);
        value_T  = unstables_for_flick(i,3);
        value_Tf = unstables_for_flick(i,4);

        [flickering(i,1), f1] = Flickering(S(:,1), value_M);
        [flickering(i,2), f2] = Flickering(S(:,2), value_Mf);
        [flickering(i,3), f3] = Flickering(S(:,3), value_T);
        [flickering(i,4), f4] = Flickering(S(:,4), value_Tf);

    if flickering(i,1)==0
        varianza_sf(i,1) = varianza(i,1);
        index{i,1} = NaN;
    else
        index{i,1} = f1;
        varianza_sf(i,1) = var(S(isfinite(S(f1,1)),1));
    end
    if flickering(i,2)==0
        varianza_sf(i,2) = varianza(i,2);
        index{i,2} = NaN;
    else
        index{i,2} = f2;
        varianza_sf(i,2) = var(S(isfinite(S(f2,2)),2));
    end
    if flickering(i,3)==0
        varianza_sf(i,3) = varianza(i,3);
        index{i,3} = NaN;
    else
        index{i,3} = f3;
        varianza_sf(i,3) = var(S(isfinite(S(f3,3)),3));
    end
    if flickering(i,4)==0
        varianza_sf(i,4) = varianza(i,4);
        index{i,4} = NaN;
    else
        index{i,4} = f4;
        varianza_sf(i,4) = var(S(isfinite(S(f4,4)),4));
    end

%%
%Graficamos las trayectorias

    if estado==1
    
        figure(2)
        subplot(3,3,10-point_ini_cond)
        hold on
        plot (T, S(:,1),'Color','k');%'#EDB120'
        title(string(10-point_ini_cond))

        line([0,nPeriods*dt], [unstables_for_flick(point_ini_cond,1), ...
            unstables_for_flick(point_ini_cond,1)], 'Color', 'b', 'LineStyle', '-')
        line([0,nPeriods*dt], [y(1),y(1)], 'Color', color, 'LineStyle', '-')
        
        xlim([0 3000])
        ylim([0, 2.5e9])
        
        figure(3)
        subplot(3,3,10-point_ini_cond)
        hold on
        plot (T, S(:,2),'Color','k');%'#EDB120'
    
        title(string(10-point_ini_cond))
        
        line([0,nPeriods*dt], [y(2),y(2)], 'Color', color, 'LineStyle', '-')
        line([0,nPeriods*dt], [unstables_for_flick(point_ini_cond,2), ...
            unstables_for_flick(point_ini_cond,2)], 'Color', 'b', 'LineStyle', '-')
        xlim([0 3000])
        ylim([0, 2.5])

        figure(4)
        subplot(3,3,10-point_ini_cond)
        hold on
        plot (T, S(:,3),'Color','k');%'#EDB120'
        title(string(10-point_ini_cond))
    
        line([0,nPeriods*dt], [y(3),y(3)], 'Color', color, 'LineStyle', '-')
        line([0,nPeriods*dt], [unstables_for_flick(point_ini_cond,3), ...
            unstables_for_flick(point_ini_cond,3)], 'Color', 'b', 'LineStyle', '-')
        xlim([0 3000])
        ylim([0,20000])

        figure(5)
        subplot(3,3,10-point_ini_cond)
        hold on
        plot (T, S(:,4),'Color','k');%'#EDB120'
        title(string(10-point_ini_cond))
    
        line([0,nPeriods*dt], [y(4),y(4)], 'Color', color, 'LineStyle', '-')
        line([0,nPeriods*dt], [unstables_for_flick(point_ini_cond,4), ...
            unstables_for_flick(point_ini_cond,4)], 'Color', 'b', 'LineStyle', '-')
        xlim([0 3000])
      %ylim([0 0.04])
%        ylim([0,20000])


else

    figure(2)
    hold on
    subplot(3,3,point_ini_cond)
    axis square
    plot (T, S(:,1),'Color','k');
    title(string(point_ini_cond))
    
    line([0,nPeriods*dt], [unstables_for_flick(point_ini_cond,1), ...
        unstables_for_flick(point_ini_cond,1)], 'Color', 'b', 'LineStyle', '-')%#4DBEEE'%#80B3FF
    line([0,nPeriods*dt], [y(1),y(1)], 'Color', color, 'LineStyle', '--')
    xlim([0 3000])
    ylim([0, 15.5e6])


        figure(3)
    hold on
    subplot(3,3,point_ini_cond)
    axis square
    plot (T, S(:,2),'Color','k');
    title(string(point_ini_cond))

    line([0,nPeriods*dt], [unstables_for_flick(point_ini_cond,2), ...
        unstables_for_flick(point_ini_cond,2)], 'Color', 'b', 'LineStyle', '-')%#4DBEEE'%#80B3FF
    line([0,nPeriods*dt], [y(2),y(2)], 'Color', color, 'LineStyle', '--')
    xlim([0 3000])
    ylim([0, 15.5e6])



    figure(4)
    hold on
    subplot(3,3,point_ini_cond)

    axis square
    plot (T, S(:,3),'Color','k');
    title(string(point_ini_cond))

    line([0,nPeriods*dt], [unstables_for_flick(point_ini_cond,3), ...
   unstables_for_flick(point_ini_cond,3)], 'Color', 'b', 'LineStyle', '-')
    line([0,nPeriods*dt], [y(3),y(3)], 'Color', color, 'LineStyle', '--')
    xlim([0 3000])
    ylim([5e6, 20.5e6])
    
    figure(5)
    hold on
    subplot(3,3,point_ini_cond)
    axis square
    plot (T, S(:,4),'Color','k');
    title(string(point_ini_cond))

    line([0,nPeriods*dt], [unstables_for_flick(point_ini_cond,4), ...
   unstables_for_flick(point_ini_cond,4)], 'Color', 'b', 'LineStyle', '-')%#4DBEEE'%#db7093
    line([0,nPeriods*dt], [y(4),y(4)], 'Color', color, 'LineStyle', '--')
    xlim([0 3000])
    ylim([-1, 10])  


    end
    
else % Si no existe estado estable
        autocor{i,1} = NaN;
        autocor{i,2} = NaN;
        autocor{i,3} = NaN;
        autocor{i,4} = NaN;

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
end % Termina filtro si no existe estado estable

clear M_sf Mf_sf T_sf Tf_sf
clear cond
end %Termina loop sobre valores de fagocitosis

%%
if estado == 1
    figure(2)
    sgtitle({'Trayectorias estocásticas','estado estable de salud'})
    han = axes(figure(2), 'Visible', 'off');

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    xlabel(han, 'tiempo');
    ylabel(han, 'M');

    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes' ...
        '\Auto_guardadas\Salud\M'];
        name='TrayectoriasM';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

    figure(3)
    sgtitle({'Trayectorias estocásticas','estado estable de salud'})
    han = axes(figure(3), 'Visible', 'off');

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    xlabel(han, 'tiempo');
    ylabel(han, 'Mf');

    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes' ...
        '\Auto_guardadas\Salud\Mf'];
        name='TrayectoriasMf';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

    figure(4)
    sgtitle({'Trayectorias estocásticas','estado estable de salud'})
    han = axes(figure(4), 'Visible', 'off');

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    xlabel(han, 'tiempo');
    ylabel(han, 'T');

    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes' ...
        '\Auto_guardadas\Salud\T'];
        name='TrayectoriasT';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

        figure(5)
    sgtitle({'Trayectorias estocásticas','estado estable de salud'})
    han = axes(figure(5), 'Visible', 'off');

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    xlabel(han, 'tiempo');
    ylabel(han, 'Tf');

    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes' ...
        '\Auto_guardadas\Salud\Tf'];
        name='TrayectoriasTf';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

else

        figure(2)
    sgtitle({'Trayectorias estocásticas','estado estable de enfermedad'})
    han = axes(figure(2), 'Visible', 'off');

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    xlabel(han, 'tiempo');
    ylabel(han, 'M');

    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes' ...
        '\Auto_guardadas\Enfermedad\M'];
        name='TrayectoriasM';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

    figure(3)
    sgtitle({'Trayectorias estocásticas','estado estable de enfermedad'})
    han = axes(figure(3), 'Visible', 'off');

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    xlabel(han, 'tiempo');
    ylabel(han, 'Mf');

    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes' ...
        '\Auto_guardadas\Enfermedad\Mf'];
        name='TrayectoriasMf';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

    figure(4)
    sgtitle({'Trayectorias estocásticas','estado estable de enfermedad'})
    han = axes(figure(4), 'Visible', 'off');

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    xlabel(han, 'tiempo');
    ylabel(han, 'T');

    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes' ...
        '\Auto_guardadas\Enfermedad\T'];
        name='TrayectoriasT';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

    figure(5)
    sgtitle({'Trayectorias estocásticas','estado estable de enfermedad'})
    han = axes(figure(5), 'Visible', 'off');

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    xlabel(han, 'tiempo');
    ylabel(han, 'Tf');

    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes' ...
        '\Auto_guardadas\Enfermedad\Tf'];
        name='TrayectoriasTf';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

end
%%
%%

au_cor = [];
for c=1:1:4
    for r = 1:1:length(autocor)
        A=autocor{r,c};
        try
        au_cor(r,c)=A(2);
        catch
        au_cor(r,c)=NaN;
        end
    end
end
%}



%% Graficámos las medidas
if estado == 1

M_autocor = rmmissing(y(:,1));
    figure
hold on
sgtitle('Autocorrelación')
subplot(4,1,1)

plot(log10(epsilon(epsilon_index)), au_cor(epsilon_index,1), '*-', 'Color', [0.5 0 0.5])
text(log10(epsilon(epsilon_index)), au_cor(epsilon_index,1), string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(au_cor(epsilon_index,1)), ...
    max(au_cor(epsilon_index,1))], 'Color', 'k', 'LineStyle', '--') ;%pto. de bif
    ylabel('autocor(M)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

subplot(4,1,2)

plot(log10(epsilon(epsilon_index)), au_cor(epsilon_index,2), '*-', 'Color', [0.5 0 0.5])
text(log10(epsilon(epsilon_index)), au_cor(epsilon_index,2), string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(au_cor(epsilon_index,2)), ...
    max(au_cor(epsilon_index,2))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('autocor(Mf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

subplot(4,1,3)

plot(log10(epsilon(epsilon_index)), au_cor(epsilon_index,3), '*-', 'Color', [0.5 0 0.5])
text(log10(epsilon(epsilon_index)), au_cor(epsilon_index,3), string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(au_cor(epsilon_index,3)), ...
    max(au_cor(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('autocor(T)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

 subplot(4,1,4)
 
plot(log10(epsilon(epsilon_index)), au_cor(epsilon_index,4), '*-', 'Color', [0.5 0 0.5])
text(log10(epsilon(epsilon_index)), au_cor(epsilon_index,4), string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(au_cor(epsilon_index,4)), ...
    max(au_cor(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    ylabel('autocor(Tf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Salud'];
        name='Autocorr';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');
    
    figure
hold on
sgtitle('Varianza')
subplot(8,1,[1 2])

plot(log10(epsilon(epsilon_index)), varianza(epsilon_index,1), '*-r')
text(log10(epsilon(epsilon_index)), varianza(epsilon_index,1), string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(epsilon_index,1)), ...
    max(varianza(epsilon_index,1))], 'Color', 'k', 'LineStyle', '--') ;%pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(M)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

subplot(8,1,[3 4])

plot(log10(epsilon(epsilon_index)), varianza(epsilon_index,2), '*-r')
text(log10(epsilon(epsilon_index)), varianza(epsilon_index,2),string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(epsilon_index,2)), ...
    max(varianza(epsilon_index,2))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(Mf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

subplot(8,1,[5 6])

plot(log10(epsilon(epsilon_index)), varianza(epsilon_index,3), '*-r')
text(log10(epsilon(epsilon_index)), varianza(epsilon_index,3),string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(epsilon_index,3)), ...
    max(varianza(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(T)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

 subplot(8,1,[7 8])
 
plot(log10(epsilon(epsilon_index)), varianza(epsilon_index,4), '*-r')
text(log10(epsilon(epsilon_index)), varianza(epsilon_index,4),string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza(epsilon_index,4)), ...
    max(varianza(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(Tf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Salud'];
        name='Varianza';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');



    figure
    hold on
sgtitle('Varianza sin flickering')
subplot(8,1,[1 2])

plot(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,1), '*-', 'Color', [0.48 0.26 0])
text(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,1),string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza_sf(epsilon_index,1)), ...
    max(varianza_sf(epsilon_index,1))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(M)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)


subplot(8,1,[3 4])
plot(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,2), '*-', 'Color', [0.48 0.26 0])
text(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,2),string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza_sf(epsilon_index,2)), ...
    max(varianza_sf(epsilon_index,2))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(Mf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)


subplot(8,1,[5 6])
plot(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,3), '*-', 'Color', [0.48 0.26 0])
text(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,3),string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza_sf(epsilon_index,3)), ...
    max(varianza_sf(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(T)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)


subplot(8,1,[7 8])
plot(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,4), '*-', 'Color', [0.48 0.26 0])
text(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,4),string(point_ini_cond:-1:1))
line([log10(bif_izquierda), log10(bif_izquierda)], [min(varianza_sf(epsilon_index,4)), ...
    max(varianza_sf(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(Tf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)
    

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Salud'];
        name='Varianza sin flickering';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

  

figure
hold on
sgtitle('Flickering')
        subplot(8,1,[1 2])
        plot(log10(epsilon(epsilon_index)), flickering(epsilon_index,1),'-g')
        hold on
        scatter(log10(epsilon(epsilon_index)), flickering(epsilon_index,1), 25,color_flick(epsilon_index,:)/255, 'filled')
        text(log10(epsilon(epsilon_index)), flickering(epsilon_index,1),string(point_ini_cond:-1:1))
        line([log10(bif_izquierda), log10(bif_izquierda)], [-1, ...
    1], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        %xlim([-8.75, -8])
        ylabel('M)', 'FontSize', 11)
        xlabel('log10(p2)', 'FontSize', 11)
    

        subplot(8,1,[3 4])
        plot(log10(epsilon(epsilon_index)), flickering(epsilon_index,2),'-g')
        hold on
        scatter(log10(epsilon(epsilon_index)), flickering(epsilon_index,2), 25,color_flick(epsilon_index,:)/255, 'filled')
        text(log10(epsilon(epsilon_index)), flickering(epsilon_index,2),string(point_ini_cond:-1:1))
       
        ylabel('Mf', 'FontSize', 11)
        xlabel('log10(p2)', 'FontSize', 11)
    

        subplot(8,1,[5 6])
        plot(log10(epsilon(epsilon_index)), flickering(epsilon_index,3),'-g')
        hold on
        scatter(log10(epsilon(epsilon_index)), flickering(epsilon_index,3), 25,color_flick(epsilon_index,:)/255, 'filled')
        text(log10(epsilon(epsilon_index)), flickering(epsilon_index,3),string(point_ini_cond:-1:1))
       
        ylabel('T', 'FontSize', 11)
        xlabel('log10(p2)', 'FontSize', 11)
    

        subplot(8,1,[7 8])
        plot(log10(epsilon(epsilon_index)), flickering(epsilon_index,4),'-g')
        hold on
        scatter(log10(epsilon(epsilon_index)), flickering(epsilon_index,4), 25,color_flick(epsilon_index,:)/255, 'filled')
        text(log10(epsilon(epsilon_index)), flickering(epsilon_index,4),string(point_ini_cond:-1:1))
     
        ylabel('Tf', 'FontSize', 11)
        xlabel('log10(p2)', 'FontSize', 11)

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Salud'];
        name='Flickering';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

figure
hold on
sgtitle('Asimetría')
    subplot(8,1,[1 2])
    plot(log10(epsilon(epsilon_index)),skew(epsilon_index,1),'*-b')
    text(log10(epsilon(epsilon_index)), skew(epsilon_index,1),string(point_ini_cond:-1:1))
    line([log10(bif_izquierda), log10(bif_izquierda)], [min(skew(epsilon_index,1)), ...
    max(skew(epsilon_index,1))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    
    ylabel('M', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)


    subplot(8,1,[3 4])
    plot(log10(epsilon(epsilon_index)),skew(epsilon_index,2),'*-b')
    text(log10(epsilon(epsilon_index)), skew(epsilon_index,2),string(point_ini_cond:-1:1))
    line([log10(bif_izquierda), log10(bif_izquierda)], [min(skew(epsilon_index,2)), ...
    max(skew(epsilon_index,2))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    
    ylabel('Mf', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)


    subplot(8,1,[5 6])
    plot(log10(epsilon(epsilon_index)),skew(epsilon_index,3),'*-b')
    text(log10(epsilon(epsilon_index)), skew(epsilon_index,3),string(point_ini_cond:-1:1))
    line([log10(bif_izquierda), log10(bif_izquierda)], [min(skew(epsilon_index,3)), ...
    max(skew(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    
    ylabel('T', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)


    subplot(8,1,[7 8])
    plot(log10(epsilon(epsilon_index)),skew(epsilon_index,4),'*-b')
    text(log10(epsilon(epsilon_index)), skew(epsilon_index,4),string(point_ini_cond:-1:1))
    line([log10(bif_izquierda), log10(bif_izquierda)], [min(skew(epsilon_index,4)), ...
    max(skew(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
   
    ylabel('T', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Salud'];
        name='Asimetría';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

else


        figure
hold on
sgtitle('Autocorrelación')

subplot(2,1,1)
plot(log10(epsilon(epsilon_index)), au_cor(epsilon_index,3), '*-', 'Color', [0.5 0 0.5])
text(log10(epsilon(epsilon_index)), au_cor(epsilon_index,3), string(1:1:point_ini_cond))
line([log10(2.98e-7), log10(2.98e-7)], [min(au_cor(epsilon_index,3)), ...
    max(au_cor(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('autocor(T)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

 subplot(2,1,2)
plot(log10(epsilon(epsilon_index)), au_cor(epsilon_index,4), '*-', 'Color', [0.5 0 0.5])
text(log10(epsilon(epsilon_index)), au_cor(epsilon_index,4), string(1:1:point_ini_cond))
line([log10(2.98e-7), log10(2.98e-7)], [min(au_cor(epsilon_index,4)), ...
    max(au_cor(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('autocor(Tf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Enfermedad'];
        name='Autocorr';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

    figure
    hold on
sgtitle('Varianza')

subplot(2,1,1)
plot(log10(epsilon(epsilon_index)), varianza(epsilon_index,3), '*-r')
text(log10(epsilon(epsilon_index)), varianza(epsilon_index,3),string(1:point_ini_cond))
line([log10(2.98e-7), log10(2.98e-7)], [min(varianza(epsilon_index,3)), ...
    max(varianza(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(T)', 'FontSize', 11)
        xlabel('log10(p2)', 'FontSize', 11)

subplot(2,1,2)
plot(log10(epsilon(epsilon_index)), varianza(epsilon_index,4), '*-r')
text(log10(epsilon(epsilon_index)), varianza(epsilon_index,4),string(1:point_ini_cond))
line([log10(2.98e-7), log10(2.98e-7)], [min(varianza(epsilon_index,4)), ...
    max(varianza(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(Tf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Enfermedad'];
        name='Varianza';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

     figure
     hold on
sgtitle('Varianza sin flickering')

subplot(2,1,1)
plot(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,3), '*-r')
text(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,3),string(1:point_ini_cond))
line([log10(2.98e-7), log10(2.98e-7)], [min(varianza_sf(epsilon_index,3)), ...
    max(varianza_sf(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(T)', 'FontSize', 11)
        xlabel('log10(p2)', 'FontSize', 11)

subplot(2,1,2)
plot(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,4), '*-r')
text(log10(epsilon(epsilon_index)), varianza_sf(epsilon_index,4),string(1:point_ini_cond))
line([log10(2.98e-7), log10(2.98e-7)], [min(varianza_sf(epsilon_index,4)), ...
    max(varianza_sf(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('var(Tf)', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Enfermedad'];
        name='Varianza sin flickering';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

    figure
    hold on
sgtitle('Flickering')

        subplot(2,1,1)
        plot(log10(epsilon(epsilon_index)), flickering(epsilon_index,3),'*-g')
        text(log10(epsilon(epsilon_index)), flickering(epsilon_index,3),string(1:point_ini_cond))
        line([log10(2.98e-7), log10(2.98e-7)], [min(flickering(epsilon_index,3)), ...
    max(flickering(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        xlim([-8.75, -6.3])
        ylabel('Número de flicks en T', 'FontSize', 11)
            xlabel('log10(p2)', 'FontSize', 11)

        subplot(2,1,2)
        plot(log10(epsilon(epsilon_index)), flickering(epsilon_index,4),'*-g')
        text(log10(epsilon(epsilon_index)), flickering(epsilon_index,4),string(1:point_ini_cond))
        line([log10(2.98e-7), log10(2.98e-7)], [min(flickering(epsilon_index,4)), ...
    max(flickering(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
        xlim([-8.75, -6.3])
        ylabel('Número de flicks en Tf', 'FontSize', 11)
        xlabel('log10(p2)', 'FontSize', 11)

        path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
            'Auto_guardadas\Enfermedad'];
        name='Flickering';
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, fullfile(path, name), 'svg');

figure
hold on
sgtitle('Asimetría')

    subplot(2,1,1)
    plot(log10(epsilon(epsilon_index)),skew(epsilon_index,3),'*-b')
    text(log10(epsilon(epsilon_index)), skew(epsilon_index,3),string(1:point_ini_cond))
    line([log10(2.98e-7), log10(2.98e-7)], [min(skew(epsilon_index,3)), ...
    max(skew(epsilon_index,3))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('Asimetría de T', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)

    subplot(2,1,2)
    plot(log10(epsilon(epsilon_index)),skew(epsilon_index,4),'*-b')
    text(log10(epsilon(epsilon_index)), skew(epsilon_index,4),string(1:point_ini_cond))
    line([log10(2.98e-7), log10(2.98e-7)], [min(skew(epsilon_index,4)), ...
    max(skew(epsilon_index,4))], 'Color', 'k', 'LineStyle', '--') ;% pto. de bif
    xlim([-8.75, -6.3])
    ylabel('Asimetría de Tf', 'FontSize', 11)
    xlabel('log10(p2)', 'FontSize', 11)
    path=['C:\Users\José Daniel\Documents\Tesis\Imágenes\' ...
    'Auto_guardadas\Enfermedad'];
    name='Asimetría';
    set(gcf, 'Renderer', 'painters');
    saveas(gcf, fullfile(path, name), 'svg');
end

toc
