% =============================================================================
% PASO 1: Carga de Datos Experimentales Completos
% =============================================================================
clear; clc;
% Configurar MATLAB en formato largo para mantener máxima precisión de bits
format long;

fprintf('--- PROCEDIMIENTO COMPLETO: SPLINE CÚBICO NATURAL (12 DECIMALES) ---\n');
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];
Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

fprintf('Puntos totales cargados en memoria: %d\n', length(f));

% =============================================================================
% PASO 2: Construcción y Procedimiento del Spline Cúbico Natural
% =============================================================================
f_espectro = linspace(min(f), max(f), 2000);

% Explicación del procedimiento en consola
fprintf('\n[PROCEDIMIENTO] planteamiento del problema:\n');
fprintf('  1. Se generan 29 intervalos entre los 30 puntos consecutivos.\n');
fprintf('  2. Cada intervalo i se modela mediante: S_i(f) = a_i*(f-f_i)^3 + b_i*(f-f_i)^2 + c_i*(f-f_i) + d_i\n');
fprintf('  3. Se resuelven las condiciones de continuidad para S(f), S''(f) y S"(f) en los nodos internos.\n');
fprintf('  4. Frontera Natural impuesta: Se obliga a que S"(f_1) = 0 y S"(f_30) = 0.\n');

% Ejecución del algoritmo analítico por trazadores cúbicos
pp = csape(f, Z, 'variational'); 

% Extraer la estructura matemática interna de la función por tramos (Piecewise Polynomial)
[break_points, coeficientes_tramos, num_tramos, orden, dimension] = unmkpp(pp);

fprintf('\n[RESULTADO] El algoritmo ha resuelto el sistema de 116 ecuaciones simultáneas.\n');
fprintf('Se han generado %d polinomios cúbicos locales independientes.\n\n', num_tramos);

% -----------------------------------------------------------------------------
% Desglose Analítico del Primer Segmento S_1(f)
% -----------------------------------------------------------------------------
% Extraemos los coeficientes del tramo 1 y redondeamos con una alta resolución de 12 decimales
coef_t1 = round(coeficientes_tramos(1,:), 12);

fprintf('--- Coeficientes del Primer Tramo S_1(f) (Dominio: 100 a 120 Hz) ---\n');
fprintf('  a_1 (Coeficiente Cúbico)     = %.12f\n', coef_t1(1));
fprintf('  b_1 (Coeficiente Cuadrático) = %.12f\n', coef_t1(2));
fprintf('  c_1 (Coeficiente Lineal)     = %.12f\n', coef_t1(3));
fprintf('  d_1 (Constante de Posición)  = %.12f\n', coef_t1(4));

fprintf('\nEcuación Matemática explícita para el Tramo 1:\n');
fprintf('  S_1(f) = (%.12f)*(f-100)^3 + (%.12f)*(f-100)^2 + (%.12f)*(f-100) + %.12f\n\n', ...
    coef_t1(1), coef_t1(2), coef_t1(3), coef_t1(4));

% =============================================================================
% PASO 3: Evaluación de la Curva Espectral Continua
% =============================================================================
% Evaluamos de forma matricial los 29 tramos sobre el vector continuo para graficar
Z_spline = round(fnval(pp, f_espectro), 12);

% =============================================================================
% PASO 4: Representación Gráfica Dedicada
% =============================================================================
fprintf('--- PASO 4: Generando y desplegando ventana gráfica ---\n');
figure('Color', [1 1 1], 'Position', [150, 150, 850, 550]);
hold on;

% 1. Graficar los 30 puntos experimentales de laboratorio (Círculos azules)
plot(f, Z, 'o', 'MarkerEdgeColor', [0.2 0.4 0.8], ...
    'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 5, 'DisplayName', 'Datos de Laboratorio');

% 2. Graficar la trayectoria del Spline Cúbico Natural por tramos (Línea verde continua)
plot(f_espectro, Z_spline, '-', 'Color', [1 0.5 0.8], 'LineWidth', 2, ...
    'DisplayName', 'Spline Cúbico Natural');

% Ajustes y configuración profesional de la visualización
title('Calibración de Sistema Portátil mediante Trazadores Cúbicos Naturales', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frecuencia f (Hz)', 'FontSize', 11);
ylabel('Magnitud de Impedancia |Z| (\Omega)', 'FontSize', 11);
xlim([0, 2900]);
ylim([125, 165]); % Límites fijos en zona física para contrastar con Runge

grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4);
legend('Location', 'best', 'FontSize', 10);
hold off;

fprintf('Proceso finalizado. La gráfica del spline continuo se muestra en pantalla.\n');

exportgraphics(gcf, 'splinecubicofinal.png')
