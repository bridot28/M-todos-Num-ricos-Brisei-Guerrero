% =============================================================================
% PASO 1: Carga de Datos y Configuración
% =============================================================================
clear; clc;
format long;

fprintf('--- ANÁLISIS DE ESTABILIDAD: SEGUNDA DERIVADA ANALÍTICA ---\n');
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];
Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

% =============================================================================
% PASO 2: Construcción del Spline y Derivación Analítica (1ra y 2da)
% =============================================================================
pp_Z = csape(f, Z, 'variational');

% Primera derivada analítica para re-localizar la raíz
pp_derivada1 = fnder(pp_Z, 1);

% Segunda derivada analítica aplicando el operador fnder nuevamente (u orden 2)
pp_derivada2 = fnder(pp_Z, 2);

% =============================================================================
% PASO 3: Localización exacta del Mínimo (Raíz de la 1ra Derivada)
% =============================================================================
f_malla_ultra = linspace(min(f), max(f), 100000);
derivada1_malla = fnval(pp_derivada1, f_malla_ultra);

idx_cruce = find(derivada1_malla(1:end-1) < 0 & derivada1_malla(2:end) > 0, 1);

f_min_exacto = f_malla_ultra(idx_cruce) - derivada1_malla(idx_cruce) * ...
    (f_malla_ultra(idx_cruce+1) - f_malla_ultra(idx_cruce)) / ...
    (derivada1_malla(idx_cruce+1) - derivada1_malla(idx_cruce));

% =============================================================================
% PASO 4: Evaluación de la Segunda Derivada en el Mínimo Encontrado
% =============================================================================
d2Z_df2_minimo = fnval(pp_derivada2, f_min_exacto);

fprintf('\n=================================================================\n');
fprintf('  RESULTADOS DE LA EVALUACIÓN DE ESTABILIDAD:\n');
fprintf('=================================================================\n');
fprintf('  Frecuencia del Mínimo (f0):  %.12f Hz\n', f_min_exacto);
fprintf('  Valor de d^2|Z|/df^2 en f0:  %.12f ohm/Hz^2\n', d2Z_df2_minimo);
fprintf('-----------------------------------------------------------------\n');

% Discusión automática del signo y la concavidad
if d2Z_df2_minimo > 0
    fprintf('  ANÁLISIS DEL SIGNO: POSITIVO (+)\n');
    fprintf('  CONCLUSIÓN: El punto crítico es un MÍNIMO LOCAL ESTABLE.\n');
    fprintf('              La curva es cóncava hacia arriba (Forma de U).\n');
elseif d2Z_df2_minimo < 0
    fprintf('  ANÁLISIS DEL SIGNO: NEGATIVO (-)\n');
    fprintf('  CONCLUSIÓN: El punto crítico es un MÁXIMO LOCAL INESTABLE.\n');
    fprintf('              La curva es cóncava hacia abajo.\n');
else
    fprintf('  CONCLUSIÓN: Segunda derivada nula. Punto de inflexión o silla.\n');
end
fprintf('=================================================================\n');

% =============================================================================
% PASO 5: Representación Gráfica de la Segunda Derivada
% =============================================================================
figure('Color', [1 1 1], 'Position', [150, 150, 850, 550]);
hold on;

% Curva continua de la segunda derivada a lo largo de la malla
f_malla_grafica = linspace(min(f), max(f), 2000);
derivada2_grafica = fnval(pp_derivada2, f_malla_grafica);
plot(f_malla_grafica, derivada2_grafica, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 2, ...
    'DisplayName', 'd^2|Z|/df^2 analítica');

% Destacar el valor exacto evaluado en el mínimo f0
plot(f_min_exacto, d2Z_df2_minimo, 'ro', 'MarkerSize', 10, 'LineWidth', 2, ...
    'MarkerFaceColor', 'r', 'DisplayName', sprintf('Evaluación en f0 (%.2f Hz)', f_min_exacto));

% Línea horizontal de referencia en Y = 0 para evaluar visualmente el signo
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'Umbral de estabilidad (Y = 0)');

% Configuración estética de los ejes
title('Segunda Derivada Analítica: Verificación de Concavidad', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frecuencia f (Hz)', 'FontSize', 11);
ylabel('d^2|Z|/df^2 (\Omega / Hz^2)', 'FontSize', 11);
xlim([0, 2900]);

grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4);
legend('Location', 'best', 'FontSize', 10);
hold off;

fprintf('\nGráfica de concavidad generada exitosamente.\n');

exportgraphics(gcf, 'segundaderivada.png')
