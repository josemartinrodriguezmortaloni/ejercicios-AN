%% EJERCICIO 1 - DERIVACIÓN NUMÉRICA
% Análisis Numérico - TP2
% Derivación de la función seno en un punto

clear all; close all; clc;

%% CONFIGURACIÓN INICIAL
% Función a analizar: y = sen(x)
% Punto de evaluación: x = π/10
% Derivada exacta: cos(π/10) ≈ 0.951056516

x0 = pi/10;  % Punto de evaluación
f = @(x) sin(x);  % Función original
df_exact = @(x) cos(x);  % Primera derivada exacta
d2f_exact = @(x) -sin(x);  % Segunda derivada exacta

% Valor exacto de la primera derivada en x0
derivada_exacta_1 = df_exact(x0);
fprintf('Valor exacto de la primera derivada en x=π/10: %.9f\n', derivada_exacta_1);

% Valor exacto de la segunda derivada en x0
derivada_exacta_2 = d2f_exact(x0);
fprintf('Valor exacto de la segunda derivada en x=π/10: %.9f\n', derivada_exacta_2);

%% PARTE A) DIFERENCIAS ADELANTADAS (DOS PUNTOS)
% Fórmula: f'(x) ≈ (f(x+h) - f(x))/h
% Error: E = -(h/2)*f''(ξ) = O(h)

fprintf('\n=== PARTE A: DIFERENCIAS ADELANTADAS ===\n');

% Diferentes valores de h (Δx)
h_values = [0.00001, 0.0001, 0.001, 0.01, 0.1];

% Inicializar vectores para almacenar resultados
deriv_adelante = zeros(size(h_values));
error_adelante = zeros(size(h_values));

fprintf('h\t\tDerivada Adelante\tError Adelante\t\tlog(h)\t\tlog(Error)\n');
fprintf('---\t\t----------------\t---------------\t\t------\t\t----------\n');

for i = 1:length(h_values)
    h = h_values(i);
    
    % Calcular derivada por diferencias adelantadas
    deriv_adelante(i) = (f(x0 + h) - f(x0)) / h;
    
    % Calcular error absoluto
    error_adelante(i) = abs(deriv_adelante(i) - derivada_exacta_1);
    
    % Mostrar resultados
    fprintf('%.5f\t%.9f\t\t%.2e\t\t%.1f\t\t%.3f\n', ...
        h, deriv_adelante(i), error_adelante(i), log10(h), log10(error_adelante(i)));
end

% Calcular pendiente de la recta log(error) vs log(h)
log_h = log10(h_values);
log_error = log10(error_adelante);
pendiente_adelante = polyfit(log_h, log_error, 1);
fprintf('\nPendiente (orden de convergencia): %.2f\n', pendiente_adelante(1));

%% PARTE B) DIFERENCIAS CENTRALES (TRES PUNTOS) - PRIMERA DERIVADA
% Fórmula: f'(x) ≈ (f(x+h) - f(x-h))/(2h)
% Error: E = -(h²/6)*f'''(ξ) = O(h²)

fprintf('\n=== PARTE B: DIFERENCIAS CENTRALES (PRIMERA DERIVADA) ===\n');

% Inicializar vectores para almacenar resultados
deriv_central_1 = zeros(size(h_values));
error_central_1 = zeros(size(h_values));

fprintf('h\t\tDerivada Central\tError Central\t\tlog(h)\t\tlog(Error)\n');
fprintf('---\t\t----------------\t---------------\t\t------\t\t----------\n');

for i = 1:length(h_values)
    h = h_values(i);
    
    % Calcular derivada por diferencias centrales
    deriv_central_1(i) = (f(x0 + h) - f(x0 - h)) / (2*h);
    
    % Calcular error absoluto
    error_central_1(i) = abs(deriv_central_1(i) - derivada_exacta_1);
    
    % Mostrar resultados
    fprintf('%.5f\t%.9f\t\t%.2e\t\t%.1f\t\t%.3f\n', ...
        h, deriv_central_1(i), error_central_1(i), log10(h), log10(error_central_1(i)));
end

% Calcular pendiente de la recta log(error) vs log(h)
log_error_central_1 = log10(error_central_1);
pendiente_central_1 = polyfit(log_h, log_error_central_1, 1);
fprintf('\nPendiente (orden de convergencia): %.2f\n', pendiente_central_1(1));

%% PARTE C) DIFERENCIAS CENTRALES (TRES PUNTOS) - SEGUNDA DERIVADA
% Fórmula: f''(x) ≈ (f(x-h) - 2f(x) + f(x+h))/h²
% Error: E = (h²/12)*f''''(ξ) = O(h²)

fprintf('\n=== PARTE C: DIFERENCIAS CENTRALES (SEGUNDA DERIVADA) ===\n');

% Valores de h para segunda derivada (según el ejercicio)
h_values_2 = [0.1, 0.01, 0.001];

% Inicializar vectores para almacenar resultados
deriv_central_2 = zeros(size(h_values_2));
error_central_2 = zeros(size(h_values_2));

fprintf('h\t\tDerivada Central\tValor Exacto\t\tError Absoluto\tlog(h)\t\tlog(Error)\n');
fprintf('---\t\t----------------\t-------------\t\t---------------\t------\t\t----------\n');

for i = 1:length(h_values_2)
    h = h_values_2(i);
    
    % Calcular segunda derivada por diferencias centrales
    deriv_central_2(i) = (f(x0 - h) - 2*f(x0) + f(x0 + h)) / (h^2);
    
    % Calcular error absoluto
    error_central_2(i) = abs(deriv_central_2(i) - derivada_exacta_2);
    
    % Mostrar resultados
    fprintf('%.3f\t\t%.9f\t\t%.9f\t\t%.2e\t\t%.1f\t\t%.3f\n', ...
        h, deriv_central_2(i), derivada_exacta_2, error_central_2(i), ...
        log10(h), log10(error_central_2(i)));
end

% Calcular pendiente de la recta log(error) vs log(h)
log_h_2 = log10(h_values_2);
log_error_central_2 = log10(error_central_2);
pendiente_central_2 = polyfit(log_h_2, log_error_central_2, 1);
fprintf('\nPendiente (orden de convergencia): %.2f\n', pendiente_central_2(1));

%% GRÁFICOS DE CONVERGENCIA
fprintf('\n=== GRÁFICOS DE CONVERGENCIA ===\n');

figure('Position', [100, 100, 1200, 800]);

% Gráfico 1: Diferencias Adelantadas
subplot(2,2,1);
loglog(h_values, error_adelante, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(h_values, h_values, 'r--', 'LineWidth', 1.5);  % Línea de referencia O(h)
xlabel('h (Δx)');
ylabel('Error Absoluto');
title('Convergencia - Diferencias Adelantadas');
legend('Error Adelantado', 'O(h)', 'Location', 'best');
grid on;

% Gráfico 2: Diferencias Centrales (Primera Derivada)
subplot(2,2,2);
loglog(h_values, error_central_1, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(h_values, h_values.^2, 'r--', 'LineWidth', 1.5);  % Línea de referencia O(h²)
xlabel('h (Δx)');
ylabel('Error Absoluto');
title('Convergencia - Diferencias Centrales (1ª Derivada)');
legend('Error Central', 'O(h²)', 'Location', 'best');
grid on;

% Gráfico 3: Diferencias Centrales (Segunda Derivada)
subplot(2,2,3);
loglog(h_values_2, error_central_2, 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(h_values_2, h_values_2.^2, 'r--', 'LineWidth', 1.5);  % Línea de referencia O(h²)
xlabel('h (Δx)');
ylabel('Error Absoluto');
title('Convergencia - Diferencias Centrales (2ª Derivada)');
legend('Error Central', 'O(h²)', 'Location', 'best');
grid on;

% Gráfico 4: Comparación de Métodos
subplot(2,2,4);
loglog(h_values, error_adelante, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(h_values, error_central_1, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
loglog(h_values_2, error_central_2, 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('h (Δx)');
ylabel('Error Absoluto');
title('Comparación de Métodos');
legend('Adelantadas O(h)', 'Centrales 1ª O(h²)', 'Centrales 2ª O(h²)', 'Location', 'best');
grid on;

sgtitle('Análisis de Convergencia - Derivación Numérica');

%% ANÁLISIS CONCEPTUAL
fprintf('\n=== ANÁLISIS CONCEPTUAL ===\n');

fprintf('\n1. GRADO DE EXACTITUD POLINÓMICA:\n');
fprintf('   - Diferencias Adelantadas (2 puntos): Exactas hasta grado 1\n');
fprintf('   - Diferencias Centrales 1ª derivada (3 puntos): Exactas hasta grado 2\n');
fprintf('   - Diferencias Centrales 2ª derivada (3 puntos): Exactas hasta grado 3\n');

fprintf('\n2. ORDEN DE CONVERGENCIAS OBSERVADAS:\n');
fprintf('   - Diferencias Adelantadas: %.2f (teórico: 1.0)\n', pendiente_adelante(1));
fprintf('   - Diferencias Centrales 1ª: %.2f (teórico: 2.0)\n', pendiente_central_1(1));
fprintf('   - Diferencias Centrales 2ª: %.2f (teórico: 2.0)\n', pendiente_central_2(1));

fprintf('\n3. COMPARACIÓN DE PRECISIÓN:\n');
fprintf('   Para la función y=sen(x):\n');
fprintf('   - En x=π/10: cos(π/10) ≈ %.6f\n', cos(pi/10));
fprintf('   - En x=π/1.8: cos(π/1.8) ≈ %.6f\n', cos(pi/1.8));
fprintf('   - La precisión depende del valor de f''(ξ) en cada punto\n');

fprintf('\n4. ERROR EN POLINOMIO DE GRADO 3:\n');
fprintf('   Si aplicamos diferencias centrales de 2ª derivada a un polinomio de grado 3:\n');
fprintf('   - El error sería C*h² donde C depende de f''''(ξ)\n');
fprintf('   - Para un polinomio de grado 3: f''''(x) = 0, por lo que el error sería 0\n');
fprintf('   - El método sería exacto para polinomios de grado ≤ 3\n');

fprintf('\n=== EJERCICIO COMPLETADO ===\n');
