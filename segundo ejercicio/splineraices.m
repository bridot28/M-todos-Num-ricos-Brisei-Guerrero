% =========================================================================
% Cálculo Analítico Exacto de Raíces usando SÓLO Spline Cúbico
% =========================================================================
clc; clear; close all;

% 1. Datos de la tabla
f = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, ...
     35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, ...
     60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, ...
     85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];

V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, ...
     1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, ...
    -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, ...
     0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

% 2. Construir el modelo Spline Cúbico
% La función 'spline' crea la estructura de datos con todos los polinomios
pp = spline(f, V);

% Extraemos los límites de los tramos (breaks) y sus coeficientes (coefs)
[breaks, coefs, tramos, orden, dim] = unmkpp(pp);

fprintf('==============================================================\n');
fprintf(' CÁLCULO DE CRUCES POR CERO (RAÍCES) SÓLO CON SPLINE CÚBICO\n');
fprintf('==============================================================\n\n');

% =========================================================================
% PRIMERA RAÍZ (Tramo entre 55.0 kHz y 57.5 kHz)
% =========================================================================
idx1 = find(breaks == 55.0); % Encontrar el índice del tramo
C1 = coefs(idx1, :);         % Extraer coeficientes [a, b, c, d]
paso = 2.5;                  % El tamaño del intervalo

% La ecuación en este tramo es: V(dx) = a*dx^3 + b*dx^2 + c*dx + d
% Donde dx = (f - 55.0)
raices_polinomio_1 = roots(C1); % Calcula las 3 raíces (pueden ser imaginarias)

% Filtrar: Buscamos la raíz real que caiga dentro de nuestro tramo [0, 2.5]
dx_real_1 = raices_polinomio_1(imag(raices_polinomio_1) == 0); 
dx_valido_1 = dx_real_1(dx_real_1 >= 0 & dx_real_1 <= paso);
raiz_1 = 55.0 + dx_valido_1; % Sumar la frecuencia base

fprintf('--- PRIMER CRUCE POR CERO ---\n');
fprintf('Tramo analizado: [%.1f, %.1f] kHz\n', breaks(idx1), breaks(idx1+1));
fprintf('Polinomio del Spline:\n');
fprintf('  V(f) = (%.6f)*dx^3 + (%.6f)*dx^2 + (%.6f)*dx + (%.6f)\n', C1(1), C1(2), C1(3), C1(4));
fprintf('  (Donde dx = f - 55.0)\n');
fprintf('Frecuencia exacta calculada: %.6f kHz\n\n', raiz_1);


% =========================================================================
% SEGUNDA RAÍZ (Tramo entre 62.5 kHz y 65.0 kHz)
% =========================================================================
idx2 = find(breaks == 62.5); % Encontrar el índice del tramo
C2 = coefs(idx2, :);         % Extraer coeficientes [a, b, c, d]

% La ecuación en este tramo es: V(dx) = a*dx^3 + b*dx^2 + c*dx + d
% Donde dx = (f - 62.5)
raices_polinomio_2 = roots(C2);

% Filtrar la raíz válida
dx_real_2 = raices_polinomio_2(imag(raices_polinomio_2) == 0);
dx_valido_2 = dx_real_2(dx_real_2 >= 0 & dx_real_2 <= paso);
raiz_2 = 62.5 + dx_valido_2; % Sumar la frecuencia base

fprintf('--- SEGUNDO CRUCE POR CERO ---\n');
fprintf('Tramo analizado: [%.1f, %.1f] kHz\n', breaks(idx2), breaks(idx2+1));
fprintf('Polinomio del Spline:\n');
fprintf('  V(f) = (%.6f)*dx^3 + (%.6f)*dx^2 + (%.6f)*dx + (%.6f)\n', C2(1), C2(2), C2(3), C2(4));
fprintf('  (Donde dx = f - 62.5)\n');
fprintf('Frecuencia exacta calculada: %.6f kHz\n\n', raiz_2);

% =========================================================================
% GRÁFICA RESUMEN
% =========================================================================
figure('Name', 'Raíces Analíticas del Spline', 'Position', [200, 200, 900, 500]);
hold on; grid on;

f_curva = linspace(50, 70, 1000);
plot(f_curva, ppval(pp, f_curva), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Spline Cúbico');
plot(f, V, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Datos Originales');
yline(0, 'r-', 'LineWidth', 1.2, 'DisplayName', 'Nivel de Alarma (V=0)');

% Marcar las raíces extraídas analíticamente
plot(raiz_1, 0, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', sprintf('Raíz 1: %.4f kHz', raiz_1));
plot(raiz_2, 0, 'g*', 'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', sprintf('Raíz 2: %.4f kHz', raiz_2));

title('Detección Exacta de Cruces por Cero (Extracción Polinomial de Spline)');
xlabel('Frecuencia f (kHz)'); ylabel('Voltaje V(f) (V)');
xlim([52 68]); ylim([-0.1 0.15]);
legend('Location', 'best');
hold off;

exportgraphics(gcf, 'splineraices.png')