% =========================================================================
% Algoritmo de Bisección sobre Spline Cúbico para AMBOS Cruces por Cero
% =========================================================================
clc; clear; close all;

% 1. Cargar todos los datos
f = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, ...
    35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, ...
    60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, ...
    85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];

V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, ...
    1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, ...
    -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, ...
    0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

% Construimos el modelo continuo mediante Spline
pp_V = spline(f, V); 

% 2. Parámetros y definición de intervalos
% Fila 1: Intervalo de la 1ra raíz | Fila 2: Intervalo de la 2da raíz
intervalos = [55.0, 57.5; 
              62.5, 65.0];

raices = zeros(1, 2); % Vector para guardar los resultados
tol = 1e-6;           % Tolerancia de error
max_iter = 30;        % Límite de seguridad

% 3. Ciclo para calcular ambas raíces
for i = 1:2
    a = intervalos(i, 1);
    b = intervalos(i, 2);
    iter = 0;
    raiz_encontrada = false;
    
    fprintf('--- INICIANDO MÉTODO DE BISECCIÓN: RAÍZ %d ---\n', i);
    fprintf('Buscando raíz en el intervalo [%.1f, %.1f] kHz\n', a, b);
    fprintf('Iter |      a      |      b      |      c (medio)  |    V(c) \n');
    fprintf('-------------------------------------------------------------------\n');

    while iter < max_iter
        iter = iter + 1;
        
        c = (a + b) / 2; % Calcular punto medio
        
        % Evaluar la función Spline
        Va = ppval(pp_V, a);
        Vc = ppval(pp_V, c);

        % Imprimir estado de la iteración
        fprintf('%2d   | %9.5f | %9.5f | %13.5f | %10.6f\n', iter, a, b, c, Vc);

        % Condición de parada
        if abs(Vc) < tol || (b - a)/2 < tol
            raiz_encontrada = true;
            break;
        end

        % Seleccionar el nuevo sub-intervalo
        if sign(Vc) == sign(Va)
            a = c;
        else
            b = c;
        end
    end

    if raiz_encontrada
        fprintf('-------------------------------------------------------------------\n');
        fprintf('¡Cruce por cero %d detectado!\n', i);
        fprintf('Frecuencia: %.4f kHz\n\n', c);
        raices(i) = c; % Guardar la raíz para la gráfica
    else
        fprintf('El método no convergió en el número máximo de iteraciones.\n\n');
    end
end

% 4. Gráfica (Zoom abarcando ambas raíces)
f_zoom = linspace(54, 66, 1000); % Rango ampliado de 54 a 66 kHz
V_zoom = ppval(pp_V, f_zoom);

figure('Name', 'Cruces por Cero - Bisección', 'Position', [200, 200, 900, 500]);
hold on; grid on;

% Dibujar la curva Spline en la zona
plot(f_zoom, V_zoom, 'b-', 'LineWidth', 2, 'DisplayName', 'Curva Spline de Voltaje');

% Dibujar los puntos originales de la tabla que caen en esa zona
puntos_f = [55.0, 57.5, 62.5, 65.0];
puntos_V = ppval(pp_V, puntos_f);
plot(puntos_f, puntos_V, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Datos de Tabla');

% Línea del cruce por cero
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Nivel Cero');

% Dibujar las raíces encontradas
plot(raices(1), 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 15, 'DisplayName', sprintf('1ra Raíz: %.4f kHz', raices(1)));
plot(raices(2), 0, 'gp', 'MarkerFaceColor', 'g', 'MarkerSize', 15, 'DisplayName', sprintf('2da Raíz: %.4f kHz', raices(2)));

title('Detección de los Cruces por Cero (Spline + Bisección)');
xlabel('Frecuencia f (kHz)');
ylabel('Voltaje V(f) (V)');
legend('Location', 'best');
xlim([54.0, 66.0]);
ylim([-0.08, 0.05]);
hold off;

exportgraphics(gcf, 'metododebiseccionraices.png')