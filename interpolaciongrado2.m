% =============================================================================
% PASO 1: Carga de Datos Experimentales Completos
% =============================================================================
clear; clc;
% Forzar a MATLAB a utilizar y mostrar su máxima precisión de punto flotante
format long;

fprintf('--- PASO 1: Cargando datos del experimento ---\n');
f_completo = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];
Z_completo = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];
fprintf('Se han cargado con éxito %d puntos de laboratorio.\n\n', length(f_completo));


% =============================================================================
% PASO 2: Selección de Nodos Clave para Grado 2 (12 decimales)
% =============================================================================
fprintf('--- PASO 2: Seleccionando los 3 nodos de soporte ---\n');
f_nodos = [100.0, 730.0, 2730.0];
Z_nodos = round([152.3, 130.9, 159.2], 12);


for i = 1:3
    fprintf('  Nodo %d: f = %.12f Hz, |Z| = %.12f ohm\n', i-1, f_nodos(i), Z_nodos(i));
end
fprintf('\n');


% =============================================================================
% PASO 3: DESARROLLO DEL MÉTODO MATRICIAL (Vandermonde)
% =============================================================================
fprintf('--- PASO 3: Procesando Método Matricial ---\n');
V = [f_nodos(:).^2, f_nodos(:), ones(3,1)];
fprintf('Matriz de Vandermonde (V):\n');
disp(V);


% Resolver el sistema lineal V * C = Z_nodos usando el operador backslash (\)
coeficientes_mat = V \ Z_nodos(:);


% Forzar redondeo estricto a 12 decimales para conservar los coeficientes pequeños
a = round(coeficientes_mat(1), 12);
b = round(coeficientes_mat(2), 12);
c = round(coeficientes_mat(3), 12);


fprintf('Coeficientes calculados (a, b, c) a 12 decimales:\n');
fprintf('  a = %.12f\n', a);
fprintf('  b = %.12f\n', b);
fprintf('  c = %.12f\n', c);
fprintf('Ecuación Matricial: P_2(f) = (%.12f)f^2 + (%.12f)f + (%.12f)\n\n', a, b, c);


% =============================================================================
% PASO 4: DESARROLLO DE LA INTERPOLACIÓN DE LAGRANGE
% =============================================================================
fprintf('--- PASO 4: Procesando Método de Lagrange ---\n');
% Cálculo y redondeo de los denominadores analíticos a 12 decimales
den0 = round((f_nodos(1) - f_nodos(2)) * (f_nodos(1) - f_nodos(3)), 12);
den1 = round((f_nodos(2) - f_nodos(1)) * (f_nodos(2) - f_nodos(3)), 12);
den2 = round((f_nodos(3) - f_nodos(1)) * (f_nodos(3) - f_nodos(2)), 12);


fprintf('  Denominador base L0: %.12f\n', den0);
fprintf('  Denominador base L1: %.12f\n', den1);
fprintf('  Denominador base L2: %.12f\n\n', den2);


% Vector de frecuencias continuo para evaluar suavemente el espectro de la gráfica
f_espectro = linspace(min(f_completo), max(f_completo), 500);


% Construcción y evaluación paso a paso del polinomio analítico de Lagrange (12 dec.)
L0 = round(((f_espectro - f_nodos(2)) .* (f_espectro - f_nodos(3))) / den0, 12);
L1 = round(((f_espectro - f_nodos(1)) .* (f_espectro - f_nodos(3))) / den1, 12);
L2 = round(((f_espectro - f_nodos(1)) .* (f_espectro - f_nodos(2))) / den2, 12);


Z_curva_lagrange = round(Z_nodos(1)*L0 + Z_nodos(2)*L1 + Z_nodos(3)*L2, 12);


% =============================================================================
% PASO 5: EVALUACIÓN DEL MÉTODO MATRICIAL Y VALIDACIÓN
% =============================================================================
fprintf('--- PASO 5: Generando y validando curvas continuas ---\n');
Z_curva_matricial = round(a * f_espectro.^2 + b * f_espectro + c, 12);


% Comparación directa del error de redondeo entre ambos caminos
error_tolerancia = max(abs(Z_curva_matricial - Z_curva_lagrange));
fprintf('Máxima discrepancia entre ambos métodos en el rango: %.12f\n', error_tolerancia);

if error_tolerancia <= 1e-11
    fprintf('VALIDACIÓN EXITOSA: Ambos métodos producen el mismo polinomio a 12 decimales.\n\n');
else
    fprintf('Aviso: Existen variaciones menores por truncamiento de punto flotante.\n\n');
end


% =============================================================================
% PASO 6: Visualización Gráfica Final en MATLAB
% =============================================================================
fprintf('--- PASO 6: Dibujando gráfico de control ---\n');
figure('Color', [1 1 1]);
hold on;


% 1. Datos reales medidos en el laboratorio (Puntos azules)
plot(f_completo, Z_completo, 'o', 'MarkerEdgeColor', [0.2 0.4 0.8], ...
     'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 5, 'DisplayName', 'Muestras Experimentales');


% 2. Nodos de control seleccionados para la interpolación (Cuadrados rojos)
plot(f_nodos, Z_nodos, 's', 'MarkerEdgeColor', 'r', ...
     'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Nodos de Control Seleccionados');


% 3. Curva continua ajustada por Lagrange/Matricial (Línea naranja)
plot(f_espectro, Z_curva_lagrange, '-', 'Color', [1 0.5 0.8], ...
     'LineWidth', 2, 'DisplayName', 'Polinomio Interpolante P_2(f) (12 dec.)');


% Formateo estético de la gráfica
title('Calibración de Bioimpedancia: Métodos de Interpolación (Grado 2)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frecuencia f (Hz)', 'FontSize', 11);
ylabel('Magnitud de Impedancia |Z| (\Omega)', 'FontSize', 11);
xlim([0, 3000]);
ylim([125, 165]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4);
legend('upper right');
hold off;


fprintf('Gráfico desplegado con éxito.\n');
exportgraphics(gcf, 'interpolaciongrado2.png')