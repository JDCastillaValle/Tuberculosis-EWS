clear 
close all
%Coeficientes
a = 1;
b = 0;
B = 1;

%Rango del parámetro
r = linspace(-1,1,400);

N_all = [];
B_all = [];

for k = 1:length(r)
    % Ecuación: -a*x^3 - b*x^2 + Bx -r = 0
    % Ecuación: -1*x^3 - 0*x^2 + 1x -r = 0
    coef = [-a -b B r(k)];
    rootsN = roots(coef);

    % Nos quedamos solo con raíces reales
    realRoots = rootsN(abs(imag(rootsN)) < 1e-6);


    x_all = [x_all; realRoots];
    B_all = [B_all; r(k)*ones(length(realRoots),1)];
end


%%
x_up = x_all(x_all >= 0.61811);
x_middle = N_all((x_all <= 0.61811) & (x_all >= -0.605964));
x_middle = sort(x_middle);
x_down = x_all(x_all <= -0.605964);
% Gráfica
figure
plot(B_all(N_all >= 0.61811), x_up, 'k',LineWidth=1.5)
hold on
plot(B_all(x_all <= -0.605964), x_down, 'k',LineWidth=1.5)
plot(B_all((x_all <= 0.61811) & (x_all >= -0.605964)), -x_middle, 'k', LineStyle='--',LineWidth=1.5)
%xlabel('p^-                                      p^+')
%ylabel('X_-                         X_c                        X_+')

text(-1.1,1,'x_+')
text(-1.1,0,'x_c')
text(-1.1,-1,'x_-')

text(-0.383459,-1.6,'p^-')
text(0.383459,-1.6,'p^+')
ylabel('p','Rotation',0,'Position',[1,-1.7]);
xlabel('x^*','Rotation',0,'Position',[-1.05,1.55]);
plot([-1,0],[0,0],'k',LineStyle='--')
plot([-1,0],[1,1],'k',LineStyle='--')
plot([-1,0],[-1,-1],'k',LineStyle='--')
plot([0.383459,0.383459],[-1.5,-0.55],'k',LineStyle="--")
plot([-0.383459,-0.383459],[-1.5,0.55],'k',LineStyle="--")
set(gca, 'XTick', [], 'YTick', [])
hold off
box off
%%
%Coeficientes
coef = [-a -b B 0];

x = linspace(-1.4, 1.4, 500);
y = polyval(coef, x); % Evaluate the polynomial
figure;
plot(x, y, 'k', 'LineWidth', 2);

hold on;

coef = [-a -b B .39];

idx = find((x>-1) & (x<0));

x_d = x(idx);
%y_d = -1*(x-2).*(x-5).*(x-8)+.53;
y_d = polyval(coef, x); % Evaluate the polynomial

plot(x_d, y_d(idx), 'k--', 'LineWidth', 1.5); % línea punteada

yline(0, 'k')


p_minus = -1;
pc = 0;
p_plus = 1;

% Marcar puntos sobre el eje
plot([p_minus pc p_plus], [0 0 0], 'ko', 'MarkerFaceColor', 'k');

% Etiquetas
text(p_minus, -0.25, 'x^-', 'HorizontalAlignment', 'center');
text(pc, -0.25, 'x^c', 'HorizontalAlignment', 'center');
text(p_plus, -0.25, 'x^+', 'HorizontalAlignment', 'center');


quiver(.6, -0.1, 0.3, 0, 0, 'k', 'MaxHeadSize', 0.5);

quiver(1.4, -0.1, -0.3, 0, 0, 'k', 'MaxHeadSize', 0.5);

quiver(-1.4, -0.1, 0.3, 0, 0, 'k', 'MaxHeadSize', 0.5);

quiver(-0.6, -0.1, -0.3, 0, 0, 'k', 'MaxHeadSize', 0.5);


ylabel('f','Rotation',0,'Position',[-1.6,1.4]);
text(1.6, 0, 'x', 'HorizontalAlignment', 'center');
text(-1.6,0,'0')

set(gca, 'XTick', [], 'YTick', [])

xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

hold off
ax = gca;
ax.Box = 'off';     % quita arriba y derecha
ax.XColor = 'none'; % oculta el eje inferior (línea y ticks)
xlabel('x','Rotation',0,'Position',[1.6,.1]);
