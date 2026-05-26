% =============================================================================
% PASO 1: Carga de Datos y Configuración
% =============================================================================
clear; clc;
format long;

fprintf('--- DERIVACIÓN ANALÍTICA MEDIANTE SPLINE CÚBICO NATURAL ---\n');
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];
Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

% =============================================================================
% PASO 2: Construcción del Spline y Derivación Analítica
% =============================================================================
% Ajuste del spline cúbico natural (variational)
pp_Z = csape(f, Z, 'variational');

% Derivación analítica del objeto spline usando la función nativa fnder
pp_derivada = fnder(pp_Z, 1);

% Calcular el valor exacto de la derivada en los 30 puntos de datos originales
derivadas_en_puntos = round(fnval(pp_derivada, f), 12);

% Mostrar los valores de la derivada en la consola para cada nodo
fprintf('[RESULTADO] Primera derivada d|Z|/df calculada en los nodos:\n');
for i = 1:length(f)
    fprintf('  f = %4d Hz  ->  d|Z|/df = %15.12f\n', f(i), derivadas_en_puntos(i));
end

% =============================================================================
% PASO 3: Localización de Precisión del Mínimo (Cruce por Cero)
% =============================================================================
% Evaluamos la derivada en una malla ultra-fina de 100,000 puntos para buscar el cero exacto
f_malla_ultra = linspace(min(f), max(f), 100000);
derivada_malla = fnval(pp_derivada, f_malla_ultra);

% Encontrar el índice donde la derivada cambia de signo (de negativa a positiva)
idx_cruce = find(derivada_malla(1:end-1) < 0 & derivada_malla(2:end) > 0, 1);

% Interpolación lineal local en la malla para hallar la raíz exacta (frecuencia del mínimo)
f_min_exacto = f_malla_ultra(idx_cruce) - derivada_malla(idx_cruce) * ...
    (f_malla_ultra(idx_cruce+1) - f_malla_ultra(idx_cruce)) / ...
    (derivada_malla(idx_cruce+1) - derivada_malla(idx_cruce));

fprintf('\n=================================================================\n');
fprintf('  UBICACIÓN DEL MÍNIMO LOCAL (Raíz de la derivada):\n');
fprintf('  Frecuencia exacta (f0): %.12f Hz\n', f_min_exacto);
fprintf('=================================================================\n');

% =============================================================================
% PASO 4: Representación Gráfica de la Derivada
% =============================================================================
figure('Color', [1 1 1], 'Position', [150, 150, 850, 550]);
hold on;

% Curva continua de la derivada analítica en la malla fina
f_malla_grafica = linspace(min(f), max(f), 2000);
derivada_grafica = fnval(pp_derivada, f_malla_grafica);
plot(f_malla_grafica, derivada_grafica, '-', 'Color', [0.0 0.4 0.7], 'LineWidth', 2, ...
     'DisplayName', 'd|Z|/df analítica');

% Puntos individuales correspondientes a las derivadas en los nodos de datos
plot(f, derivadas_en_puntos, 'o', 'MarkerEdgeColor', [0.2 0.5 0.8], ...
     'MarkerFaceColor', [1 1 1], 'MarkerSize', 5, 'DisplayName', 'Derivadas en Nodos');

% Línea de referencia horizontal en Y = 0 (Eje de cambio de signo)
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'Cruce por Cero (d|Z|/df = 0)');

% Destacar la ubicación exacta del mínimo experimental
plot(f_min_exacto, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2, ...
     'DisplayName', sprintf('Mínimo: %.4f Hz', f_min_exacto));

% Configuración del gráfico
title('Primera Derivada Analítica de la Impedancia con Spline Cúbico', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frecuencia f (Hz)', 'FontSize', 11);
ylabel('d|Z|/df (\Omega / Hz)', 'FontSize', 11);
xlim([0, 2900]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4);
legend('Location', 'best', 'FontSize', 10);
hold off;

exportgraphics(gcf, 'primeraderivada.png')
