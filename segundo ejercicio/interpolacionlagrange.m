% Script de MATLAB para Interpolación de Lagrange de 2do Grado
clc; clear; close all;

% --- CASO 1 y 2: Estimaciones para f = 41.0 kHz ---
% Los 3 puntos más cercanos son f = 37.5, 40.0, 42.5
f_41 = [37.5, 40.0, 42.5];
V_41 = [1.048, 0.866, 0.689];
Z_41 = [145.2, 145.8, 147.3];
x_est_41 = 41.0;

% --- CASO 3 y 4: Estimaciones para f = 73.0 kHz ---
% Los 3 puntos más cercanos son f = 70.0, 72.5, 75.0
f_73 = [70.0, 72.5, 75.0];
V_73 = [0.197, 0.318, 0.452];
Z_73 = [200.1, 203.1, 205.2];
x_est_73 = 73.0;

% Crear figura para las gráficas
figure('Name', 'Interpolación de Lagrange', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 800]);

% 1. Calcular V(41.0 kHz)
subplot(2,2,1);
calcular_lagrange(f_41, V_41, x_est_41, 'Voltaje V(f) a 41 kHz');

% 2. Calcular |Z|(41.0 kHz)
subplot(2,2,2);
calcular_lagrange(f_41, Z_41, x_est_41, 'Impedancia |Z|(f) a 41 kHz');

% 3. Calcular V(73.0 kHz)
subplot(2,2,3);
calcular_lagrange(f_73, V_73, x_est_73, 'Voltaje V(f) a 73 kHz');

% 4. Calcular |Z|(73.0 kHz)
subplot(2,2,4);
calcular_lagrange(f_73, Z_73, x_est_73, 'Impedancia |Z|(f) a 73 kHz');

% =========================================================================
% FUNCIÓN PARA CALCULAR, MOSTRAR PASOS Y GRAFICAR LAGRANGE
% =========================================================================
function calcular_lagrange(x_pts, y_pts, x_est, titulo)
syms x;

% 1. Cálculo de los Polinomios base de Lagrange (L0, L1, L2)
% Fórmula: L_i(x) = Productoria de (x - x_j) / (x_i - x_j)
L0 = (x - x_pts(2))*(x - x_pts(3)) / ((x_pts(1) - x_pts(2))*(x_pts(1) - x_pts(3)));
L1 = (x - x_pts(1))*(x - x_pts(3)) / ((x_pts(2) - x_pts(1))*(x_pts(2) - x_pts(3)));
L2 = (x - x_pts(1))*(x - x_pts(2)) / ((x_pts(3) - x_pts(1))*(x_pts(3) - x_pts(2)));

% 2. Polinomio de interpolación de 2do grado final
% Fórmula: P(x) = y0*L0 + y1*L1 + y2*L2
P = y_pts(1)*L0 + y_pts(2)*L1 + y_pts(3)*L2;

% Expandir y simplificar el polinomio para obtener la ecuación de la forma ax^2 + bx + c
% Utilizamos vpa para limitar a 6 decimales significativos y que sea legible
P_simp = vpa(expand(P), 6);

% 3. Evaluar el polinomio en el punto solicitado
val_est = double(subs(P, x, x_est));

% --- IMPRESIÓN EN CONSOLA (Paso a Paso) ---
fprintf('======================================================\n');
fprintf('Análisis para: %s\n', titulo);
fprintf('Puntos utilizados (x, y):\n');
fprintf(' P0 = (%.1f, %.3f)\n', x_pts(1), y_pts(1));
fprintf(' P1 = (%.1f, %.3f)\n', x_pts(2), y_pts(2));
fprintf(' P2 = (%.1f, %.3f)\n', x_pts(3), y_pts(3));
fprintf('\nEcuación del polinomio de Lagrange simplificada:\n P(x) = %s\n', char(P_simp));
fprintf('\nCálculo de la estimación en x = %.1f:\n Valor estimado = %.4f\n', x_est, val_est);
fprintf('======================================================\n\n');

% --- CREACIÓN DE LA GRÁFICA ---
hold on; grid on;
% Rango de la gráfica (un poco más amplio que los puntos dados)
rango_x = [min(x_pts)-1, max(x_pts)+1];

% Dibujar la curva del polinomio
fplot(P, rango_x, 'b-', 'LineWidth', 1.5);

% Dibujar los 3 puntos base que construyen el polinomio
plot(x_pts, y_pts, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% Dibujar el punto calculado/estimado
plot(x_est, val_est, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

title(titulo);
xlabel('Frecuencia f (kHz)');
ylabel('Valor (V o \Omega)');
legend('Polinomio P(x)', 'Datos tabla', 'Valor estimado', 'Location', 'best');
hold off;
exportgraphics(gcf, 'interpolacionlagrange.png')
end