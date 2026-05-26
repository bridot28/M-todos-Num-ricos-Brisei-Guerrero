% =========================================================================
% Spline Cúbico Natural para toda la tabla de mediciones
% =========================================================================
clc; clear; close all;

% 1. Extraer todos los datos de la tabla (40 mediciones provistas)
f_full = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, ...
          35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, ...
          60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, ...
          85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];

V_full = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, ...
          1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, ...
         -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, ...
          0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

Z_full = [182.4, 178.9, 175.1, 171.0, 166.8, 162.7, 158.9, 155.4, 152.0, 149.0, ...
          146.1, 145.2, 145.8, 147.3, 149.9, 153.5, 158.0, 163.2, 168.9, 174.8, ...
          180.5, 186.2, 191.5, 196.2, 200.1, 203.1, 205.2, 206.3, 206.1, 204.7, ...
          198.0, 194.4, 190.9, 187.8, 185.1, 183.0, 181.6, 180.8, 180.6, 180.9];

% 2. Calcular los coeficientes de los Splines para V y |Z|
[a_v, b_v, c_v, d_v] = coeffs_spline_natural(f_full, V_full);
[a_z, b_z, c_z, d_z] = coeffs_spline_natural(f_full, Z_full);

% 3. Imprimir ecuaciones de los tramos específicos y calcular estimaciones
disp('CÁLCULOS CON SPLINE CÚBICO NATURAL:');
V_41 = analizar_y_mostrar(f_full, a_v, b_v, c_v, d_v, 41.0, 'Voltaje V(41.0 kHz)');
Z_41 = analizar_y_mostrar(f_full, a_z, b_z, c_z, d_z, 41.0, 'Impedancia |Z|(41.0 kHz)');
V_73 = analizar_y_mostrar(f_full, a_v, b_v, c_v, d_v, 73.0, 'Voltaje V(73.0 kHz)');
Z_73 = analizar_y_mostrar(f_full, a_z, b_z, c_z, d_z, 73.0, 'Impedancia |Z|(73.0 kHz)');

% 4. Graficar los resultados de toda la curva
f_curva = linspace(min(f_full), max(f_full), 1000);
V_curva = eval_spline(f_full, a_v, b_v, c_v, d_v, f_curva);
Z_curva = eval_spline(f_full, a_z, b_z, c_z, d_z, f_curva);

figure('Name', 'Spline Cúbico Natural Completo', 'NumberTitle', 'off', 'Position', [100, 100, 900, 700]);

% Gráfica de Voltaje
subplot(2,1,1); hold on; grid on;
plot(f_full, V_full, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Datos Tabla');
plot(f_curva, V_curva, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Spline Cúbico (Toda la curva)');
plot([41.0, 73.0], [V_41, V_73], 'r*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Estimaciones Solicitadas');
title('Spline Cúbico Natural - Voltaje V(f)');
xlabel('Frecuencia f (kHz)'); ylabel('Voltaje (V)');
legend('Location', 'best');

% Gráfica de Impedancia
subplot(2,1,2); hold on; grid on;
plot(f_full, Z_full, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Datos Tabla');
plot(f_curva, Z_curva, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Spline Cúbico (Toda la curva)');
plot([41.0, 73.0], [Z_41, Z_73], 'r*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Estimaciones Solicitadas');
title('Spline Cúbico Natural - Impedancia |Z|(f)');
xlabel('Frecuencia f (kHz)'); ylabel('Impedancia (\Omega)');
legend('Location', 'best');

% =========================================================================
% FUNCIONES AUXILIARES PARA EL CÁLCULO DEL SPLINE MATRICIAL
% =========================================================================

% Función para calcular los coeficientes [a, b, c, d] de cada tramo
function [a, b, c, d] = coeffs_spline_natural(x, y)
    n = length(x); h = diff(x);
    alpha = zeros(1, n-1);
    for i = 2:n-1
        alpha(i) = 3/h(i)*(y(i+1)-y(i)) - 3/h(i-1)*(y(i)-y(i-1));
    end
    l = zeros(1, n); mu = zeros(1, n); z = zeros(1, n); l(1) = 1;
    for i = 2:n-1
        l(i) = 2*(x(i+1)-x(i-1)) - h(i-1)*mu(i-1);
        mu(i) = h(i)/l(i);
        z(i) = (alpha(i)-h(i-1)*z(i-1))/l(i);
    end
    l(n) = 1; c_temp = zeros(1, n); b = zeros(1, n-1); d = zeros(1, n-1); a = y(1:n-1);
    for j = n-1:-1:1
        c_temp(j) = z(j) - mu(j)*c_temp(j+1);
        b(j) = (y(j+1)-y(j))/h(j) - h(j)*(c_temp(j+1)+2*c_temp(j))/3;
        d(j) = (c_temp(j+1)-c_temp(j))/(3*h(j));
    end
    c = c_temp(1:n-1);
end

% Función para evaluar la curva dadas las ecuaciones
function yq = eval_spline(x, a, b, c, d, xq)
    yq = zeros(size(xq));
    n = length(x);
    for k = 1:length(xq)
        idx = find(x <= xq(k), 1, 'last');
        if isempty(idx), idx = 1; end
        if idx >= n, idx = n-1; end
        dx = xq(k) - x(idx);
        yq(k) = a(idx) + b(idx)*dx + c(idx)*dx^2 + d(idx)*dx^3;
    end
end

% Función para imprimir el análisis de un punto y mostrar su polinomio
function val = analizar_y_mostrar(f, a, b, c, d, f_est, titulo)
    idx = find(f <= f_est, 1, 'last');
    if idx == length(f); idx = length(f)-1; end
    
    dx = f_est - f(idx);
    val = a(idx) + b(idx)*dx + c(idx)*dx^2 + d(idx)*dx^3;
    
    % Mostrar ecuación simplificada usando variables simbólicas
    syms x_sym
    S_eq = a(idx) + b(idx)*(x_sym - f(idx)) + c(idx)*(x_sym - f(idx))^2 + d(idx)*(x_sym - f(idx))^3;
    S_eq_expanded = vpa(expand(S_eq), 5);
    
    fprintf('------------------------------------------------------\n');
    fprintf('%s\n', titulo);
    fprintf('El punto %.1f kHz cae en el tramo interpolado [%.1f, %.1f]\n', f_est, f(idx), f(idx+1));
    fprintf('Ecuación en este tramo:\n S(f) = %s\n', char(S_eq_expanded));
    fprintf('Valor estimado final = %.4f\n', val);
    fprintf('------------------------------------------------------\n');
exportgraphics(gcf, 'splinegraficos.png')
end