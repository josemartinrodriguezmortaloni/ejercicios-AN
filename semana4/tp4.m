% TP Runge-Kutta 2do Orden - Script completo

% Función principal para resolver todos los ejercicios
function tp_runge_kutta()
    fprintf('=== TP RUNGE-KUTTA 2DO ORDEN ===\n\n');
    
    % Ejercicio 2
    ejercicio2();
    
    % Ejercicio 3
    ejercicio3();
    
    % Ejercicio 4
    ejercicio4();
    
    % Ejercicio 7 (sistema)
    ejercicio7();
end

% Método RK2 general
function [t, y] = rk2(f, tspan, y0, dt, omega)
    t = tspan(1):dt:tspan(2);
    n = length(t);
    m = length(y0);
    y = zeros(m, n);
    y(:,1) = y0;
    
    for i = 1:n-1
        k1 = dt * f(t(i), y(:,i));
        yG = y(:,i) + k1/(2*omega);
        tG = t(i) + dt/(2*omega);
        k2 = dt * f(tG, yG);
        y(:,i+1) = y(:,i) + (1-omega)*k1 + omega*k2;
    end
end

% Ejercicio 2
function ejercicio2()
    fprintf('EJERCICIO 2:\n');
    f = @(t,u) 2*u - 2*t - 1;
    u_exact = @(t) exp(2*t) + t + 1;
    
    dt = 0.25;
    [t_mej, u_mej] = rk2(f, [0 0.5], 2, dt, 0.5);  % Euler Mejorado
    [t_mod, u_mod] = rk2(f, [0 0.5], 2, dt, 1.0);  % Euler Modificado
    
    fprintf('Euler Mejorado: u(0.5) = %.6f\n', u_mej(end));
    fprintf('Euler Modificado: u(0.5) = %.6f\n', u_mod(end));
    fprintf('Solución exacta: u(0.5) = %.6f\n', u_exact(0.5));
    fprintf('Error Mejorado: %.6f\n', abs(u_mej(end) - u_exact(0.5)));
    fprintf('Error Modificado: %.6f\n\n', abs(u_mod(end) - u_exact(0.5)));
    
    % Gráfica con dt fino
    [t_graf, u_graf] = rk2(f, [0 10], 2, 0.01, 0.5);
    t_ex = 0:0.01:10;
    
    figure(1);
    subplot(2,1,1);
    plot(t_graf, u_graf, 'b-', t_ex, u_exact(t_ex), 'r--');
    xlabel('t'); ylabel('u(t)');
    title('Ejercicio 2: Solución de EDO');
    legend('RK2', 'Exacta');
    grid on;
    
    subplot(2,1,2);
    plot(t_graf, abs(u_graf - u_exact(t_graf)), 'g-');
    xlabel('t'); ylabel('Error absoluto');
    title('Error del método');
    grid on;
end

% Ejercicio 3
function ejercicio3()
    fprintf('EJERCICIO 3:\n');
    f = @(t,u) -u/2 + t/2;
    u_exact = @(t) 6*exp(-t/2) - 2 + t;
    
    dt = 0.25;
    [t_mej, u_mej] = rk2(f, [0 0.5], 4, dt, 0.5);
    [t_mod, u_mod] = rk2(f, [0 0.5], 4, dt, 1.0);
    
    % Método de Euler simple para comparar
    [t_euler, u_euler] = euler_simple(f, [0 0.5], 4, dt);
    
    fprintf('RK2 Mejorado: u(0.5) = %.6f\n', u_mej(end));
    fprintf('RK2 Modificado: u(0.5) = %.6f\n', u_mod(end));
    fprintf('Euler simple: u(0.5) = %.6f\n', u_euler(end));
    fprintf('Exacta: u(0.5) = %.6f\n\n', u_exact(0.5));
    
    % Gráficas comparativas
    figure(2);
    dt_fino = 0.01;
    [t_fino, u_fino] = rk2(f, [0 10], 4, dt_fino, 0.5);
    t_ex = 0:0.01:10;
    
    plot(t_fino, u_fino, 'b-', t_ex, u_exact(t_ex), 'r--');
    xlabel('t'); ylabel('u(t)');
    title('Ejercicio 3: Comparación de soluciones');
    legend('RK2', 'Exacta');
    grid on;
end

% Ejercicio 4
function ejercicio4()
    fprintf('EJERCICIO 4:\n');
    f = @(t,y) -2*y + exp(-2*t);
    y_exact = @(t) 0.1*exp(-2*t) + t*exp(-2*t);
    
    dt = 0.2;
    [t, y] = rk2(f, [0 0.4], 0.1, dt, 0.5);
    
    fprintf('y(0.4) aproximado = %.6f\n', y(end));
    fprintf('y(0.4) exacto = %.6f\n', y_exact(0.4));
    fprintf('Error = %.6f\n\n', abs(y(end) - y_exact(0.4)));
    
    % Buscar dt para error < 1%
    y_max = max(abs(y_exact(0:0.01:6)));
    for dt_test = [0.1 0.05 0.02 0.01 0.005]
        [t_test, y_test] = rk2(f, [0 6], 0.1, dt_test, 0.5);
        t_exact_test = 0:dt_test:6;
        error_norm = max(abs(y_test - y_exact(t_exact_test)));
        fprintf('dt = %.3f -> Error norma inf = %.3f%% del máximo\n', ...
                dt_test, 100*error_norm/y_max);
    end
    fprintf('\n');
end

% Ejercicio 7 (Sistema)
function ejercicio7()
    fprintf('EJERCICIO 7 (Sistema):\n');
    A = [-10 4; -4 0];
    f = @(t, x) A*x;
    x0 = [5; 3];
    
    dt = 0.01;
    [t, x] = rk2(f, [0 0.02], x0, dt, 0.5);
    
    fprintf('Con dt=%.3f:\n', dt);
    fprintf('x1(0.01) = %.4f, x2(0.01) = %.4f\n', x(1,2), x(2,2));
    fprintf('x1(0.02) = %.4f, x2(0.02) = %.4f\n\n', x(1,3), x(2,3));
    
    % Solución para tiempo mayor
    [t_largo, x_largo] = rk2(f, [0 2], x0, 0.02, 0.5);
    
    figure(3);
    subplot(1,2,1);
    plot(t_largo, x_largo(1,:), 'b-', t_largo, x_largo(2,:), 'r-');
    xlabel('t'); ylabel('x(t)');
    title('Evolución temporal');
    legend('x1(t)', 'x2(t)');
    grid on;
    
    subplot(1,2,2);
    plot(x_largo(1,:), x_largo(2,:), 'g-');
    xlabel('x1'); ylabel('x2');
    title('Diagrama de fase');
    grid on;
end

% Método de Euler simple (para comparación)
function [t, y] = euler_simple(f, tspan, y0, dt)
    t = tspan(1):dt:tspan(2);
    n = length(t);
    y = zeros(size(y0,1), n);
    y(:,1) = y0;
    
    for i = 1:n-1
        y(:,i+1) = y(:,i) + dt * f(t(i), y(:,i));
    end
end

% Ejecutar el TP
tp_runge_kutta();