% =============================================================================
% PASO 1: Carga de Datos y Configuración del Umbral
% =============================================================================
clear; clc;
format long;

fprintf('--- CARACTERIZACIÓN DE BANDA SEGURA: RAÍCES DE |Z|(f) = 150 ohm ---\n\n');
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];
Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

Z_th = 150.0; % Umbral crítico de atenuación

% =============================================================================
% PASO 2: Construcción del Spline y Desplazamiento del Umbral
% =============================================================================
pp = csape(f, Z, 'variational');

% Desarmar el polinomio a tramos
[break_points, coeficientes_tramos, num_tramos, orden, dimension] = unmkpp(pp);

raices_encontradas = [];

% =============================================================================
% PASO 3: Búsqueda Analítica de Raíces Tramo por Tramo
% =============================================================================
for i = 1:num_tramos
    % Límites del tramo actual
    f_izq = break_points(i);
    f_der = break_points(i+1);
    
    % Coeficientes originales del tramo: [a, b, c, d]
    coef_originales = coeficientes_tramos(i, :);
    
    % Para resolver S_i(f) - 150 = 0, le restamos Z_th al término constante (d)
    coef_desplazados = coef_originales;
    coef_desplazados(4) = coef_desplazados(4) - Z_th;
    
    % roots() calcula las raíces de la ecuación: a*dx^3 + b*dx^2 + c*dx + (d - 150) = 0
    % donde dx = f - f_izq
    raices_locales_dx = roots(coef_desplazados);
    
    % Evaluar cada una de las raíces encontradas en este tramo
    for r = 1:length(raices_locales_dx)
        % Filtro 1: Nos quedamos solo con las raíces reales (descartar complejas conjugadas)
        if isreal(raices_locales_dx(r))
            % Convertir la variable local dx de vuelta a frecuencia absoluta f
            f_raiz = f_izq + raices_locales_dx(r);
            
            % Filtro 2: Verificar que la raíz caiga estrictamente dentro del tramo i
            % Usamos una pequeña tolerancia por redondeo numérico (1e-9)
            if (f_raiz >= f_izq - 1e-9) && (f_raiz <= f_der + 1e-9)
                % Almacenar la raíz válida encontrada
                raices_encontradas = [raices_encontradas, f_raiz];
            end
        end
    end
end

% Eliminar posibles raíces duplicadas en las fronteras exactas de los nodos
raices_encontradas = unique(round(raices_encontradas, 9));

% =============================================================================
% PASO 4: Impresión de Resultados en Consola
% =============================================================================
fprintf('=================================================================\n');
fprintf('        FRECUENCIAS LÍMITE DE LA BANDA DE OPERACIÓN SEGURA\n');
fprintf('=================================================================\n');
for idx = 1:length(raices_encontradas)
    fprintf('  Frecuencia Límite %d (f_raíz): %.12f Hz\n', idx, raices_encontradas(idx));
end
fprintf('-----------------------------------------------------------------\n');
if length(raices_encontradas) >= 2
    f_inf = raices_encontradas(1);
    f_sup = raices_encontradas(2);
    fprintf('  ANÁLISIS DE LA BANDA SEGURA (|Z| < 150 ohm):\n');
    fprintf('  Rango de operación libre de atenuación: [%.4f Hz — %.4f Hz]\n', f_inf, f_sup);
    fprintf('  Ancho de banda seguro efectivo (BW):    %.4f Hz\n', f_sup - f_inf);
else
    fprintf('  Aviso: No se completaron las dos raíces esperadas en el rango.\n');
end
fprintf('=================================================================\n');

% =============================================================================
% PASO 5: Visualización Gráfica del Umbral de Atenuación
% =============================================================================
f_espectro = linspace(min(f), max(f), 2000);
Z_spline = fnval(pp, f_espectro);

figure('Color', [1 1 1], 'Position', [150, 150, 850, 550]);
hold on;

% 1. Curva del Spline Cúbico de Impedancia
plot(f_espectro, Z_spline, '-', 'Color', [0.0 0.5 0.3], 'LineWidth', 2, 'DisplayName', 'Spline Cúbico |Z|(f)');

% 2. Línea horizontal del Umbral Crítico Z_th = 150
yline(Z_th, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Umbral de Atenuación (Z_{th} = 150 \Omega)');

% 3. Destacar las raíces encontradas con marcadores de estrella negros
for idx = 1:length(raices_encontradas)
    plot(raices_encontradas(idx), Z_th, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
         'DisplayName', sprintf('Corte %d: %.2f Hz', idx, raices_encontradas(idx)));
end

% Ajustes estéticos
title('Identificación de Frecuencias Límite para la Banda de Operación', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frecuencia f (Hz)', 'FontSize', 11);
ylabel('Magnitud de Impedancia |Z| (\Omega)', 'FontSize', 11);
xlim([0, 2900]);
ylim([125, 165]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4);
legend('Location', 'best');
hold off;

% Malla local ultra-fina alrededor de la raíz superior
f_raiz2 = 2216.591398246371;
h_paso = 0.001; % Delta f

% Evaluación del spline en los entornos de la raíz
Z_mas  = fnval(pp, f_raiz2 + h_paso);
Z_menos = fnval(pp, f_raiz2 - h_paso);

% Derivada numérica de la impedancia (Diferencia Central)
dZ_df_num = (Z_mas - Z_menos) / (2 * h_paso);

% Sensibilidad (Derivada inversa numérica)
df_dZ_num = 1 / dZ_df_num;

fprintf('Sensibilidad Numérica Calculada: %.12f Hz/ohm\n', df_dZ_num);

exportgraphics(gcf, 'raices.png')
