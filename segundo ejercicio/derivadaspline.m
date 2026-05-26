% Cálculo de Derivada Analítica usando Spline Cúbico
clc; clear;

% Datos (los 40 puntos)
f = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];
V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, 1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, 0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

% Construir estructura de Spline Cúbico
pp_V = csape(f, V, 'variational'); % 'variational' = Spline natural (segunda derivada cero en bordes)

% Extraer los coeficientes polinomiales [d, c, b, a] para cada tramo
% pp_V.coefs tiene tamaño (N-1) x 4
coefs = pp_V.coefs;

% Puntos a evaluar
f_eval = [10.0, 40.0, 70.0, 100.0];

disp('--- DERIVADA CON SPLINE CÚBICO (Sensibilidad V/kHz) ---');
for i = 1:length(f_eval)
    % Encontrar en qué tramo del spline cae la frecuencia
    tramo_idx = find(f == f_eval(i));
    
    % Si estamos en el último punto, usamos el último tramo
    if tramo_idx == length(f)
        tramo_idx = length(f) - 1;
        dx = f(end) - f(end-1);
    else
        dx = 0; % Evaluando en el nodo exacto izquierdo
    end
    
    % Ecuación de la derivada: S'(x) = 3*a*dx^2 + 2*b*dx + c
    % En la matriz coefs, las columnas son [a, b, c, d]
    a = coefs(tramo_idx, 1);
    b = coefs(tramo_idx, 2);
    c = coefs(tramo_idx, 3);
    
    derivada = 3*a*dx^2 + 2*b*dx + c;
    
    fprintf('f = %.1f kHz -> dV/df = %.4f\n', f_eval(i), derivada);
end