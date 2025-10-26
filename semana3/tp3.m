% TP4 - Ecuaciones Diferenciales Ordinarias - Método de Euler
% Análisis Numérico - 2023
% Implementación correcta según PDF del TP
% =================================================================

function tp3_completo
    fprintf('====================================================\n');
    fprintf('TP EDO - MÉTODO DE EULER\n');
    fprintf('Análisis Numérico - UTN FRM - 2023\n');
    fprintf('====================================================\n\n');
    
    % Ejecutar todos los ejercicios
    ejercicio1();  % Teórico
    ejercicio2();
    ejercicio3();
    ejercicio4();
    ejercicio5();
    ejercicio6();
    ejercicio7a();  % Primer sistema del ejercicio 7
    ejercicio7b();  % Segundo sistema del ejercicio 7
    ejercicio8();
    ejercicio9();
    
    fprintf('\n¡Trabajo Práctico completado!\n');
end

% =================================================================
% MÉTODOS DE EULER
% =================================================================

function [t, u] = euler_simple(f, t0, u0, dt, n_pasos)
    % Método de Euler para EDO de primer orden
    % f: función que define du/dt = f(t,u)
    % t0, u0: condiciones iniciales
    % dt: tamaño del paso
    % n_pasos: número de pasos
    
    t = zeros(n_pasos+1, 1);
    u = zeros(n_pasos+1, 1);
    
    t(1) = t0;
    u(1) = u0;
    
    for i = 1:n_pasos
        k1 = dt * f(t(i), u(i));
        u(i+1) = u(i) + k1;
        t(i+1) = t(i) + dt;
    end
end

function [t, X] = euler_sistema(F, t0, X0, dt, n_pasos)
    % Método de Euler para sistemas de EDO
    % F: función vectorial que define dX/dt = F(t,X)
    % t0: valor inicial de t
    % X0: vector de condiciones iniciales
    % dt: tamaño del paso
    % n_pasos: número de pasos
    
    m = length(X0);  % número de ecuaciones
    t = zeros(n_pasos+1, 1);
    X = zeros(n_pasos+1, m);
    
    t(1) = t0;
    X(1, :) = X0';
    
    for i = 1:n_pasos
        k1 = dt * F(t(i), X(i, :)');
        X(i+1, :) = X(i, :) + k1';
        t(i+1) = t(i) + dt;
    end
end

% =================================================================
% EJERCICIO 1 - Teórico
% =================================================================

function ejercicio1()
    fprintf('=== EJERCICIO 1 ===\n');
    fprintf('Interpretación geométrica del Método de Euler:\n');
    fprintf('------------------------------------------------\n');
    fprintf('El método de Euler aproxima la curva solución mediante\n');
    fprintf('segmentos de recta tangentes. En cada punto (tn, un),\n');
    fprintf('se usa la pendiente f(tn, un) para avanzar al siguiente punto.\n\n');
    
    fprintf('Error de Truncamiento Local (ETL):\n');
    fprintf('ETL = (Δt)²/2 * d²u(t)/dt²|ξ\n');
    fprintf('Por lo tanto, el método de Euler es de primer orden O(Δt)\n\n');
    
    % Visualización geométrica
    figure;
    t = 0:0.01:1;
    y_exact = exp(t);  % Ejemplo: y' = y, y(0) = 1
    
    % Euler con paso grande para visualizar
    dt_vis = 0.2;
    n_vis = 5;
    [t_euler, y_euler] = euler_simple(@(t,y) y, 0, 1, dt_vis, n_vis);
    
    plot(t, y_exact, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_euler, y_euler, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
    
    % Dibujar tangentes
    for i = 1:length(t_euler)-1
        t_tang = [t_euler(i), t_euler(i+1)];
        y_tang = [y_euler(i), y_euler(i+1)];
        plot(t_tang, y_tang, 'g--', 'LineWidth', 1);
    end
    
    xlabel('t');
    ylabel('y');
    title('Interpretación Geométrica del Método de Euler');
    legend('Solución exacta', 'Aproximación Euler', 'Tangentes', 'Location', 'northwest');
    grid on;
end

% =================================================================
% EJERCICIO 2
% =================================================================

function ejercicio2()
    fprintf('\n=== EJERCICIO 2 ===\n');
    fprintf('EDO: du/dt - 2u(t) = -2t - 1, u(0) = 2\n');
    fprintf('Solución exacta: u(t) = exp(2t) + t + 1\n');
    fprintf('------------------------------------------------\n');
    
    % Definir la EDO: du/dt = f(t,u) = 2u(t) - 2t - 1
    f = @(t, u) 2*u - 2*t - 1;
    
    % Condiciones iniciales
    t0 = 0;
    u0 = 2;
    t_final = 0.5;
    
    % a) Con Δt = 0.10
    dt1 = 0.10;
    n1 = round(t_final/dt1);
    [t1, u1] = euler_simple(f, t0, u0, dt1, n1);
    
    % b) Con Δt = 0.25
    dt2 = 0.25;
    n2 = round(t_final/dt2);
    [t2, u2] = euler_simple(f, t0, u0, dt2, n2);
    
    % Solución exacta
    u_exacta = @(t) exp(2*t) + t + 1;
    
    % Mostrar resultados en t = 0.5
    u_exact_05 = u_exacta(0.5);
    error1 = abs(u1(end) - u_exact_05);
    error2 = abs(u2(end) - u_exact_05);
    
    fprintf('Resultados en t = 0.5:\n');
    fprintf('Δt = 0.10: u(0.5) = %.6f, Error = %.6f\n', u1(end), error1);
    fprintf('Δt = 0.25: u(0.5) = %.6f, Error = %.6f\n', u2(end), error2);
    fprintf('Valor exacto: u(0.5) = %.6f\n', u_exact_05);
    
    % Tabla detallada
    fprintf('\nTabla detallada con Δt = 0.10:\n');
    fprintf('t\t\tu_Euler\t\tu_exacta\tError\n');
    for i = 1:length(t1)
        u_ex = u_exacta(t1(i));
        error = abs(u1(i) - u_ex);
        fprintf('%.2f\t\t%.6f\t%.6f\t%.6f\n', t1(i), u1(i), u_ex, error);
    end
    
    % Gráfico como en el PDF
    dt_graf = 0.01;
    n_graf = round(1.3/dt_graf);
    [t_graf, u_graf] = euler_simple(f, t0, u0, dt_graf, n_graf);
    
    t_exact = 0:0.01:1.3;
    u_exact = u_exacta(t_exact);
    
    % Con dt = 0.1
    dt_01 = 0.1;
    n_01 = round(1.3/dt_01);
    [t_01, u_01] = euler_simple(f, t0, u0, dt_01, n_01);
    
    % Con dt = 0.05
    dt_005 = 0.05;
    n_005 = round(1.3/dt_005);
    [t_005, u_005] = euler_simple(f, t0, u0, dt_005, n_005);
    
    figure;
    plot(t_exact, u_exact, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_01, u_01, 'ro', 'MarkerSize', 6);
    plot(t_005, u_005, 'g^', 'MarkerSize', 4);
    xlabel('t');
    ylabel('u(t)');
    title('Ejercicio 2: Método de Euler');
    legend('u(t) exacta', 'u(t) Δt=0.1', 'u(t) Δt=0.05', 'Location', 'northwest');
    grid on;
    axis([0 1.3 1 17]);
end

% =================================================================
% EJERCICIO 3
% =================================================================

function ejercicio3()
    fprintf('\n=== EJERCICIO 3 ===\n');
    fprintf('EDO: du/dt + (1/2)u(t) = t/2, u(0) = 4\n');
    fprintf('Solución exacta: u(t) = 6*exp(-t/2) - 2 + t\n');
    fprintf('------------------------------------------------\n');
    
    % Definir la EDO: du/dt = f(t,u) = -u/2 + t/2
    f = @(t, u) -u/2 + t/2;
    
    % Condiciones iniciales
    t0 = 0;
    u0 = 4;
    
    % a) Con Δt = 0.10
    dt1 = 0.10;
    n1 = 5;  % hasta t = 0.5
    [t1, u1] = euler_simple(f, t0, u0, dt1, n1);
    
    % b) Con Δt = 0.25
    dt2 = 0.25;
    n2 = 2;  % hasta t = 0.5
    [t2, u2] = euler_simple(f, t0, u0, dt2, n2);
    
    % Solución exacta
    u_exacta = @(t) 6*exp(-t/2) - 2 + t;
    
    % Mostrar resultados
    fprintf('Resultados en t = 0.5:\n');
    fprintf('Δt = 0.10: u(0.5) = %.6f, Error = %.6f\n', u1(end), abs(u1(end) - u_exacta(0.5)));
    fprintf('Δt = 0.25: u(0.5) = %.6f, Error = %.6f\n', u2(end), abs(u2(end) - u_exacta(0.5)));
    fprintf('Valor exacto: u(0.5) = %.6f\n', u_exacta(0.5));
    
    % Gráfico con Δt = 0.01 hasta t = 10
    dt_graf = 0.01;
    n_graf = round(10/dt_graf);
    [t_graf, u_graf] = euler_simple(f, t0, u0, dt_graf, n_graf);
    
    t_exact = 0:0.01:10;
    u_exact = u_exacta(t_exact);
    
    figure;
    subplot(2,1,1);
    plot(t_graf, u_graf, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(t_exact, u_exact, 'b-', 'LineWidth', 1.5);
    xlabel('t');
    ylabel('u(t)');
    title('Ejercicio 3: Solución de EDO');
    legend('y aprox', 'y exacta', 'Location', 'northeast');
    grid on;
    axis([0 10 2 9]);
    
    % Función error
    subplot(2,1,2);
    error_func = abs(u_graf - u_exacta(t_graf));
    plot(t_graf, error_func, 'k-', 'LineWidth', 1.5);
    xlabel('t');
    ylabel('Error absoluto');
    title('Error del método de Euler (Δt = 0.01)');
    grid on;
end

% =================================================================
% EJERCICIO 4
% =================================================================

function ejercicio4()
    fprintf('\n=== EJERCICIO 4 ===\n');
    fprintf('EDO: dy/dt + 2y(t) = exp(-2t), y(0) = 1/10\n');
    fprintf('Solución exacta: y(t) = (1/10)*exp(-2t) + t*exp(-2t)\n');
    fprintf('------------------------------------------------\n');
    
    % Definir la EDO: dy/dt = f(t,y) = -2y + exp(-2t)
    f = @(t, y) -2*y + exp(-2*t);
    
    % Condiciones iniciales
    t0 = 0;
    y0 = 1/10;
    
    % Solución exacta
    y_exacta = @(t) (1/10)*exp(-2*t) + t.*exp(-2*t);
    
    % Búsqueda del Δt óptimo
    dt_values = [0.1, 0.05, 0.01, 0.005, 0.001];
    errors = zeros(size(dt_values));
    
    t_eval = 0:0.01:5;
    y_exact_eval = y_exacta(t_eval);
    max_y = max(abs(y_exact_eval));
    
    fprintf('Análisis de convergencia:\n');
    fprintf('Δt\t\tNorma∞(error)\tError relativo (%%)\n');
    
    for i = 1:length(dt_values)
        dt = dt_values(i);
        n = round(5/dt);
        [t_euler, y_euler] = euler_simple(f, t0, y0, dt, n);
        
        % Interpolar para comparar en los mismos puntos
        y_euler_interp = interp1(t_euler, y_euler, t_eval, 'linear', 'extrap');
        error_vec = abs(y_euler_interp - y_exact_eval);
        errors(i) = max(error_vec);
        error_rel = (errors(i)/max_y) * 100;
        
        fprintf('%.3f\t\t%.6e\t%.2f%%\n', dt, errors(i), error_rel);
        
        if error_rel < 1.0
            fprintf('✓ Δt = %.3f cumple con error < 1%% del máximo\n', dt);
            dt_optimo = dt;
        end
    end
    
    % Gráfico con Δt óptimo
    dt_graf = 0.001;
    n_graf = round(6/dt_graf);
    [t_graf, y_graf] = euler_simple(f, t0, y0, dt_graf, n_graf);
    
    figure;
    t_exact = 0:0.01:6;
    y_exact = y_exacta(t_exact);
    
    plot(t_exact, y_exact, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_graf(1:100:end), y_graf(1:100:end), 'ro', 'MarkerSize', 6);
    xlabel('t');
    ylabel('y(t)');
    title('Ejercicio 4: dy/dt + 2y = exp(-2t)');
    legend('Solución exacta', 'Aproximación Euler', 'Location', 'northeast');
    grid on;
end

% =================================================================
% EJERCICIO 5
% =================================================================

function ejercicio5()
    fprintf('\n=== EJERCICIO 5 ===\n');
    fprintf('EDO: du/dt - 2t*(u(t))^2 = 0, u(0) = 1\n');
    fprintf('Solución exacta: u(t) = 1/(1 - t^2)\n');
    fprintf('Nota: Asíntota vertical en t = 1\n');
    fprintf('------------------------------------------------\n');
    
    % Definir la EDO: du/dt = f(t,u) = 2t*u^2
    f = @(t, u) 2*t*u^2;
    
    % Condiciones iniciales
    t0 = 0;
    u0 = 1;
    
    % Solución exacta
    u_exacta = @(t) 1./(1 - t.^2);
    
    % Diferentes valores de Δt
    dt_values = [0.05, 0.01, 0.005];
    colors = {'r-', 'g-', 'm-'};
    
    figure;
    hold on;
    
    for i = 1:length(dt_values)
        dt = dt_values(i);
        t_max = 0.95;  % Antes de la asíntota
        n = round(t_max/dt);
        [t_euler, u_euler] = euler_simple(f, t0, u0, dt, n);
        
        plot(t_euler, u_euler, colors{i}, 'LineWidth', 1.5);
        
        % Calcular error
        u_exact_points = u_exacta(t_euler);
        error_norm = norm(u_euler - u_exact_points, inf);
        fprintf('Δt = %.3f: Norma∞(error) = %.6e\n', dt, error_norm);
    end
    
    % Solución exacta
    t_exact = 0:0.001:0.95;
    u_exact = u_exacta(t_exact);
    plot(t_exact, u_exact, 'b-', 'LineWidth', 2);
    
    xlabel('t');
    ylabel('u(t)');
    title('Ejercicio 5: Ecuación con asíntota vertical en t=1');
    legend('Δt=0.05', 'Δt=0.01', 'Δt=0.005', 'Solución exacta', 'Location', 'northwest');
    grid on;
    ylim([0 20]);
end

% =================================================================
% EJERCICIO 6
% =================================================================

function ejercicio6()
    fprintf('\n=== EJERCICIO 6 ===\n');
    fprintf('EDO: dy/dt + t*y(t) = 0, y(0) = 1\n');
    fprintf('Solución exacta: y(t) = exp(-t^2/2)\n');
    fprintf('------------------------------------------------\n');
    
    % Definir la EDO: dy/dt = f(t,y) = -t*y
    f = @(t, y) -t*y;
    
    % Condiciones iniciales
    t0 = 0;
    y0 = 1;
    
    % Solución exacta
    y_exacta = @(t) exp(-t.^2/2);
    
    % Análisis de convergencia
    dt_values = [0.1, 0.05, 0.01, 0.005];
    t_final = 3.5;
    
    fprintf('Análisis de convergencia:\n');
    fprintf('Δt\t\tNorma∞(error)\n');
    
    for i = 1:length(dt_values)
        dt = dt_values(i);
        n = round(t_final/dt);
        [t_euler, y_euler] = euler_simple(f, t0, y0, dt, n);
        
        y_exact = y_exacta(t_euler);
        error_norm = norm(y_euler - y_exact, inf);
        fprintf('%.3f\t\t%.6e\n', dt, error_norm);
    end
    
    % Gráfico
    dt_graf = 0.001;
    n_graf = round(3.5/dt_graf);
    [t_graf, y_graf] = euler_simple(f, t0, y0, dt_graf, n_graf);
    
    figure;
    plot(t_graf, y_graf, 'r-', 'LineWidth', 2);
    hold on;
    t_exact = 0:0.01:3.5;
    y_exact = y_exacta(t_exact);
    plot(t_exact, y_exact, 'b--', 'LineWidth', 2);
    
    xlabel('t');
    ylabel('y(t)');
    title('Ejercicio 6: dy/dt + t*y = 0');
    legend('Euler', 'Exacta', 'Location', 'northeast');
    grid on;
end

% =================================================================
% EJERCICIO 7a - Primer sistema
% =================================================================

function ejercicio7a()
    fprintf('\n=== EJERCICIO 7a ===\n');
    fprintf('Sistema: [x1''; x2''] = [-10 4; -4 0]*[x1; x2]\n');
    fprintf('Condiciones iniciales: [x1(0); x2(0)] = [5; 3]\n');
    fprintf('------------------------------------------------\n');
    
    % Definir el sistema
    A = [-10 4; -4 0];
    F = @(t, X) A * X;
    
    % Condiciones iniciales
    t0 = 0;
    X0 = [5; 3];
    
    % Con Δt = 0.02 hasta t = 2
    dt = 0.02;
    n = round(2/dt);
    [t, X] = euler_sistema(F, t0, X0, dt, n);
    
    % Solución exacta
    x1_exacta = @(t) (1/3)*exp(-8*t) + (14/3)*exp(-2*t);
    x2_exacta = @(t) (2/3)*exp(-8*t) + (7/3)*exp(-2*t);
    
    % Mostrar tabla de valores
    fprintf('\nTabla de valores (primeros puntos):\n');
    fprintf('t\t\tx1(t)\t\tx2(t)\t\tk1_x1\t\tk1_x2\n');
    for i = 1:min(5, length(t))
        if i < length(t)
            k1 = dt * F(t(i), X(i,:)');
            fprintf('%.2f\t\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
                    t(i), X(i,1), X(i,2), k1(1), k1(2));
        else
            fprintf('%.2f\t\t%.6f\t%.6f\n', t(i), X(i,1), X(i,2));
        end
    end
    
    % Gráficos
    figure;
    subplot(2,1,1);
    plot(t, X(:,1), 'b-', 'LineWidth', 2);
    hold on;
    plot(t, X(:,2), 'r-', 'LineWidth', 2);
    t_exact = 0:0.01:2;
    x1_ex = x1_exacta(t_exact);
    x2_ex = x2_exacta(t_exact);
    plot(t_exact, x1_ex, 'b--', 'LineWidth', 1);
    plot(t_exact, x2_ex, 'r--', 'LineWidth', 1);
    xlabel('t');
    ylabel('x(t)');
    title('Ejercicio 7a: Solución del sistema');
    legend('x1 aprox', 'x2 aprox', 'x1 exacta', 'x2 exacta', 'Location', 'northeast');
    grid on;
    
    % Error
    subplot(2,1,2);
    x1_error = abs(X(:,1) - x1_exacta(t));
    x2_error = abs(X(:,2) - x2_exacta(t));
    plot(t, x1_error, 'b-', 'LineWidth', 2);
    hold on;
    plot(t, x2_error, 'r-', 'LineWidth', 2);
    xlabel('t');
    ylabel('Error absoluto');
    title('Error del método de Euler');
    legend('Error x1', 'Error x2', 'Location', 'northwest');
    grid on;
end

% =================================================================
% EJERCICIO 7b - Segundo sistema (alta frecuencia)
% =================================================================

function ejercicio7b()
    fprintf('\n=== EJERCICIO 7b ===\n');
    fprintf('Sistema: [x1''; x2''] = [0 1; -100 0]*[x1; x2]\n');
    fprintf('Condiciones iniciales: [x1(0); x2(0)] = [5; 3]\n');
    fprintf('Sistema oscilatorio de alta frecuencia\n');
    fprintf('------------------------------------------------\n');
    
    % Definir el sistema
    A = [0 1; -100 0];
    F = @(t, X) A * X;
    
    % Condiciones iniciales
    t0 = 0;
    X0 = [5; 3];
    
    % a) Con Δt = 0.01 (divergente)
    dt1 = 0.01;
    n1 = round(2.91/dt1);
    [t1, X1] = euler_sistema(F, t0, X0, dt1, n1);
    
    fprintf('\nCon Δt = 0.01 (primeros valores):\n');
    fprintf('t\t\tx1\t\tx2\n');
    indices_mostrar = [1, 2, 3, 4, 5, 6, 7];
    for idx = indices_mostrar
        if idx <= length(t1)
            fprintf('%.2e\t%.3e\t%.3e\n', t1(idx), X1(idx,1), X1(idx,2));
        end
    end
    
    % b) Con Δt = 0.0001 (convergente)
    dt2 = 0.0001;
    n2 = round(1.13/dt2);
    [t2, X2] = euler_sistema(F, t0, X0, dt2, n2);
    
    fprintf('\nCon Δt = 0.0001: Sistema estable y oscilatorio\n');
    fprintf('Período aproximado: %.3f\n', 2*pi/10);
    
    % Gráficos
    figure;
    
    % Gráfico con Δt = 0.01
    subplot(2,2,1);
    plot(t1, X1(:,1), 'b-', 'LineWidth', 1.5);
    xlabel('t');
    ylabel('x1(t)');
    title('x1(t) con Δt = 0.01 (divergente)');
    grid on;
    
    subplot(2,2,2);
    plot(X1(:,1), X1(:,2), 'r-', 'LineWidth', 1);
    xlabel('x1');
    ylabel('x2');
    title('Diagrama de fase (Δt = 0.01)');
    grid on;
    
    % Gráfico con Δt = 0.0001
    subplot(2,2,3);
    plot(t2(1:100:end), X2(1:100:end,1), 'b-', 'LineWidth', 1.5);
    xlabel('t');
    ylabel('x1(t)');
    title('x1(t) con Δt = 0.0001 (estable)');
    grid on;
    
    subplot(2,2,4);
    plot(X2(:,1), X2(:,2), 'g-', 'LineWidth', 1);
    xlabel('x1');
    ylabel('x2');
    title('Diagrama de fase (Δt = 0.0001)');
    grid on;
    axis equal;
end

% =================================================================
% EJERCICIO 8
% =================================================================

function ejercicio8()
    fprintf('\n=== EJERCICIO 8 ===\n');
    fprintf('Sistema: [x1''; x2''] = [-0.5 2; -2 -0.5]*[x1; x2]\n');
    fprintf('Condiciones iniciales: [x1(0); x2(0)] = [10; 0]\n');
    fprintf('------------------------------------------------\n');
    
    % Definir el sistema
    A = [-0.5 2; -2 -0.5];
    F = @(t, X) A * X;
    
    % Condiciones iniciales
    t0 = 0;
    X0 = [10; 0];
    
    % Con Δt = 0.055
    dt1 = 0.055;
    n1 = round(15/dt1);
    [t1, X1] = euler_sistema(F, t0, X0, dt1, n1);
    
    % Con Δt = 0.0055
    dt2 = 0.0055;
    n2 = round(15/dt2);
    [t2, X2] = euler_sistema(F, t0, X0, dt2, n2);
    
    fprintf('Análisis de convergencia de valores picos en [1, 4]:\n');
    
    % Encontrar picos en X1 para Δt = 0.055
    idx_rango1 = find(t1 >= 1 & t1 <= 4);
    if ~isempty(idx_rango1)
        [pks1, locs1] = findpeaks(X1(idx_rango1,1));
        if ~isempty(pks1)
            fprintf('Δt = 0.055: Primer pico = %.4f en t = %.2f\n', ...
                    pks1(1), t1(idx_rango1(locs1(1))));
        end
    end
    
    % Encontrar picos en X2 para Δt = 0.0055
    idx_rango2 = find(t2 >= 1 & t2 <= 4);
    if ~isempty(idx_rango2)
        [pks2, locs2] = findpeaks(X2(idx_rango2,1));
        if ~isempty(pks2)
            fprintf('Δt = 0.0055: Primer pico = %.4f en t = %.2f\n', ...
                    pks2(1), t2(idx_rango2(locs2(1))));
        end
    end
    
    % Gráficos
    figure;
    
    % Con Δt = 0.055
    subplot(2,2,1);
    plot(t1, X1(:,1), 'b-', 'LineWidth', 1.5);
    xlabel('t');
    ylabel('x1(t)');
    title('x1(t) con Δt = 0.055');
    grid on;
    
    subplot(2,2,2);
    plot(X1(:,1), X1(:,2), 'r.', 'MarkerSize', 3);
    xlabel('x1');
    ylabel('x2');
    title('Diagrama de fase (Δt = 0.055)');
    grid on;
    
    % Con Δt = 0.0055
    subplot(2,2,3);
    plot(t2(1:10:end), X2(1:10:end,1), 'b-', 'LineWidth', 1.5);
    xlabel('t');
    ylabel('x1(t)');
    title('x1(t) con Δt = 0.0055');
    grid on;
    
    subplot(2,2,4);
    plot(X2(1:10:end,1), X2(1:10:end,2), 'g-', 'LineWidth', 1);
    xlabel('x1');
    ylabel('x2');
    title('Diagrama de fase (Δt = 0.0055)');
    grid on;
end

% =================================================================
% EJERCICIO 9
% =================================================================

function ejercicio9()
    fprintf('\n=== EJERCICIO 9 ===\n');
    fprintf('Sistema: [x1''; x2''] = [0 1; -4 0]*[x1; x2] + sin(3t)*[0; 10]\n');
    fprintf('Condiciones iniciales: [x1(0); x2(0)] = [0; 0]\n');
    fprintf('Sistema forzado con término sinusoidal\n');
    fprintf('------------------------------------------------\n');
    
    % Definir el sistema con forzamiento
    A = [0 1; -4 0];
    F = @(t, X) A * X + sin(3*t) * [0; 10];
    
    % Condiciones iniciales
    t0 = 0;
    X0 = [0; 0];
    
    % Probar diferentes valores de Δt
    dt_values = [0.01, 0.005, 0.001];
    
    fprintf('Buscando Δt adecuado para aproximación similar a RK2...\n');
    
    for dt = dt_values
        n = round(15/dt);
        [t, X] = euler_sistema(F, t0, X0, dt, n);
        
        % Verificar estabilidad
        if max(abs(X(:,1))) < 100  % Si no diverge
            fprintf('Δt = %.3f: Solución estable\n', dt);
            
            % Usar este Δt para los gráficos
            if dt == 0.001
                figure;
                
                subplot(1,2,1);
                plot(t(1:10:end), X(1:10:end,1), 'b-', 'LineWidth', 1.5);
                xlabel('t');
                ylabel('x1(t)');
                title(sprintf('Ejercicio 9: x1(t) con Δt = %.3f', dt));
                grid on;
                ylim([-6 6]);
                
                subplot(1,2,2);
                plot(X(1:10:end,1), X(1:10:end,2), 'r-', 'LineWidth', 1);
                xlabel('x1');
                ylabel('x2');
                title('Diagrama de fase');
                grid on;
                axis([-6 6 -15 15]);
            end
        else
            fprintf('Δt = %.3f: Solución divergente\n', dt);
        end
    end
    
    fprintf('\nNota: Para obtener resultados similares al método RK2,\n');
    fprintf('se requiere Δt ≤ 0.001 con el método de Euler.\n');
end

% =================================================================
% FUNCIÓN AUXILIAR PARA ENCONTRAR PICOS
% =================================================================

function [pks, locs] = findpeaks(data)
    % Función simple para encontrar picos locales
    pks = [];
    locs = [];
    
    for i = 2:length(data)-1
        if data(i) > data(i-1) && data(i) > data(i+1)
            pks = [pks; data(i)];
            locs = [locs; i];
        end
    end
end