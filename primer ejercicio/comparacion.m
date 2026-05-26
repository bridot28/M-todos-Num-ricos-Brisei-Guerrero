% =============================================================================
% PASO 1: Carga de Datos y Configuración de la Malla Fina
% =============================================================================
clear; clc;
format long;

fprintf('--- COMPARATIVA: POLINOMIO DE GRADO 2 VS. SPLINE CÚBICO NATURAL ---\n');
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];
Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

% Definición de la malla fina de frecuencia (2000 puntos densos)
f_malla_fina = linspace(min(f), max(f), 2000);

% =============================================================================
% PASO 2: Evaluación del Polinomio Seleccionado (Grado 2 - 3 Nodos)
% =============================================================================
fprintf('\n[PROCEDIMIENTO 1] Calculando Polinomio de Grado 2...\n');
f_nodos = [100.0, 730.0, 2730.0];
Z_nodos = [152.3, 130.9, 159.2];

% Matriz de Vandermonde local para los 3 nodos
V2 = [f_nodos(:).^2, f_nodos(:), ones(3,1)];
coef_p2 = V2 \ Z_nodos(:);
coef_p2 = round(coef_p2, 12); % Alta precisión

fprintf('  Ecuación P_2(f) = (%.12f)*f^2 + (%.12f)*f + (%.12f)\n', ...
    coef_p2(1), coef_p2(2), coef_p2(3));

% Evaluación en la malla fina
Z_malla_p2 = polyval(coef_p2, f_malla_fina);

% =============================================================================
% PASO 3: Evaluación del Spline Cúbico Natural (30 Puntos)
% =============================================================================
fprintf('\n[PROCEDIMIENTO 2] Calculando Spline Cúbico Natural por tramos...\n');
pp = csape(f, Z, 'variational'); 

% Evaluación en la malla fina
Z_malla_spline = fnval(pp, f_malla_fina);
fprintf('  Spline evaluado exitosamente sobre la malla fina de 2000 puntos.\n');

% =============================================================================
% PASO 4: Representación Gráfica Comparativa
% =============================================================================
fprintf('\n--- PASO 4: Generando gráfica comparativa ---\n');
figure('Color', [1 1 1], 'Position', [150, 150, 900, 600]);
hold on;

% 1. Puntos experimentales completos de laboratorio
plot(f, Z, 'o', 'MarkerEdgeColor', [0.2 0.4 0.8], ...
    'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 5, 'DisplayName', 'Datos de Laboratorio');

% 2. Nodos específicos utilizados para el polinomio cuadrático
plot(f_nodos, Z_nodos, 's', 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', '3 Nodos de Control (P_2)');

% 3. Curva de la malla fina para el Polinomio de Grado 2
plot(f_malla_fina, Z_malla_p2, '--', 'Color', [0.8 0.4 0.0], 'LineWidth', 2, ...
    'DisplayName', 'Polinomio de Grado 2 (Ajuste Global)');

% 4. Curva de la malla fina para el Spline Cúbico Natural
plot(f_malla_fina, Z_malla_spline, '-', 'Color', [1 0.5 0.8], 'LineWidth', 2, ...
    'DisplayName', 'Spline Cúbico Natural (Ajuste Local)');

% Configuración de etiquetas y diseño profesional
title('Evaluación en Malla Fina: Polinomio Global vs. Spline Local', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frecuencia f (Hz)', 'FontSize', 11);
ylabel('Magnitud de Impedancia |Z| (\Omega)', 'FontSize', 11);
xlim([0, 2900]);
ylim([125, 165]);

grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4);
legend('Location', 'best', 'FontSize', 10);
hold off;

fprintf('Procesamiento e informe gráfico finalizados.\n');

exportgraphics(gcf, 'comparacion.png')
