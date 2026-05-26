% Script para Cálculo de Derivadas Numéricas (Sensibilidad)
clc; clear; close all;

% 1. Definición de Datos
f = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, ...
    35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, ...
    60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, ...
    85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];

V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, ...
    1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, ...
    -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, ...
    0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

h = 2.5; % Tamaño de paso constante (kHz)

% 2. Puntos a evaluar
f_eval = [40.0, 70.0, 100.0];
dV_cd2 = zeros(1, 3);
dV_cd4 = zeros(1, 3);
idx_eval = zeros(1, 3);

fprintf('--- RESULTADOS DE DERIVADA dV/df (Sensibilidad) ---\n');

% 3. Cálculos de Diferencias Finitas
for i = 1:3
    % Encontrar el índice en el arreglo
    idx = find(f == f_eval(i));
    idx_eval(i) = idx;

    % Diferencia Centrada de Orden 2
    dV_cd2(i) = (V(idx+1) - V(idx-1)) / (2*h);

    % Diferencia Centrada de Orden 4
    dV_cd4(i) = (-V(idx+2) + 8*V(idx+1) - 8*V(idx-1) + V(idx-2)) / (12*h);

    % Imprimir resultados
    fprintf('\nFrecuencia: %.1f kHz\n', f_eval(i));
    fprintf('  Orden 2 [O(h^2)]: %.6f V/kHz\n', dV_cd2(i));
    fprintf('  Orden 4 [O(h^4)]: %.6f V/kHz\n', dV_cd4(i));
end

% 4. Generación de Gráfica (Visualizando la Derivada como Tangente)
figure('Name', 'Sensibilidad dV/df', 'Position', [150, 150, 800, 500]);
hold on; grid on;

% Dibujar curva base
plot(f, V, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'Color', [0.6 0.6 0.6], 'DisplayName', 'Datos V(f)');

% Colores para cada tangente
colores = {'r', 'b', 'g'};

for i = 1:3
    x0 = f_eval(i);
    y0 = V(idx_eval(i));
    m = dV_cd4(i); % Usamos la de orden 4 por ser más precisa

    % Crear línea tangente (y - y0 = m*(x - x0)) -> y = m*(x - x0) + y0
    x_tan = linspace(x0 - 5, x0 + 5, 2); % Línea que abarca +/- 5 kHz
    y_tan = m * (x_tan - x0) + y0;

    % Dibujar Tangente
    plot(x_tan, y_tan, colores{i}, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Tangente O(h^4) en %.1f kHz (m=%.4f)', x0, m));

    % Resaltar el punto
    plot(x0, y0, [colores{i} '*'], 'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

title('Estimación de la Sensibilidad dV/df (Derivada Numérica)');
xlabel('Frecuencia f (kHz)');
ylabel('Voltaje V(f) (V)');
legend('Location', 'best');
hold off;
exportgraphics(gcf, 'Primeraderivadaestimacion.png')