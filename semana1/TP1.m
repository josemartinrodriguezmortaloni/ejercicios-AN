% Trabajo Práctico N°1 - Integración Numérica
% Ejercicio 1: Regla del Trapecio

function TP1
    clear; clc; close all;

    % Función a integrar: sin(πt)
    f = @(t) sin(pi*t);

    a = 0;
    b = 0.5;

    % Valores de n (número de subdivisiones)
    valores_n = [10, 20, 40, 80];

    fprintf('=== EJERCICIO 1: REGLA DEL TRAPECIO ===\n');
    fprintf('Integración de sen(πt) en 1/4 del período\n\n');

    % ∫[0,0.5] sen(πt)dt = [-cos(πt)/π] desde 0 hasta 0.5 = 1/π
    integral_exacta = 1/pi;

    fprintf('Función: sin(πt)\n');
    fprintf('Intervalo: [%.1f, %.1f]\n', a, b);
    fprintf('Integral exacta: 1/π = %.8f\n\n', integral_exacta);

    % TABLA 1: SUMA DE RIEMANN (Im)
    fprintf('TABLA 1: SUMA DE RIEMANN (Im)\n');
    fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', 'N', 'Δt', 'Im(Δt)', 'Error', 'Log(Δt)', 'Log(Er)', 'Pendiente');
    fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', '---', '-----', '-------', '-----', '-------', '-------', '--------');

    % Vectores para almacenar resultados de Riemann
    delta_t_r = zeros(1, length(valores_n));
    errores_r = zeros(1, length(valores_n));

    for i = 1:length(valores_n)
        n = valores_n(i);
        dt = (b-a)/n;  % Δt
        resultado_riemann = suma_riemann(f, a, b, n);
        error_riemann = abs(resultado_riemann - integral_exacta);

        delta_t_r(i) = dt;
        errores_r(i) = error_riemann;

        % Calcular logaritmos
        log_dt = log10(dt);
        log_er = log10(error_riemann);

        % Calcular pendiente (orden de error) con el punto anterior
        if i > 1
            pendiente = (log_er - log10(errores_r(i-1))) / (log_dt - log10(delta_t_r(i-1)));
            fprintf('%-5d %-8.5f %-15.8f %-12.6f %-10.5f %-10.5f %-10.6f\n', n, dt, resultado_riemann, error_riemann, log_dt, log_er, pendiente);
        else
            fprintf('%-5d %-8.5f %-15.8f %-12.6f %-10.5f %-10.5f %-10s\n', n, dt, resultado_riemann, error_riemann, log_dt, log_er, '---');
        end
    end

    fprintf('\n');

    % TABLA 2: TRAPECIOS MÚLTIPLE (ITrap)
    fprintf('TABLA 2: TRAPECIOS MÚLTIPLE (ITrap)\n');
    fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', 'N', 'Δt', 'ITrap(Δt)', 'Error', 'Log(Δt)', 'Log(Er)', 'Pendiente');
    fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', '---', '-----', '-------', '-----', '-------', '-------', '--------');

    % Vectores para almacenar resultados de Trapecios
    delta_t = zeros(1, length(valores_n));
    errores = zeros(1, length(valores_n));

    for i = 1:length(valores_n)
        n = valores_n(i);
        dt = (b-a)/n;  % Δt
        resultado = regla_trapecio(f, a, b, n);
        error = abs(resultado - integral_exacta);

        delta_t(i) = dt;
        errores(i) = error;

        % Calcular logaritmos
        log_dt = log10(dt);
        log_er = log10(error);

        % Calcular pendiente (orden de error) con el punto anterior
        if i > 1
            pendiente = (log_er - log10(errores(i-1))) / (log_dt - log10(delta_t(i-1)));
            fprintf('%-5d %-8.5f %-15.8f %-12.6f %-10.5f %-10.5f %-10.6f\n', n, dt, resultado, error, log_dt, log_er, pendiente);
        else
            fprintf('%-5d %-8.5f %-15.8f %-12.6f %-10.5f %-10.5f %-10s\n', n, dt, resultado, error, log_dt, log_er, '---');
        end
    end
    fprintf('\n\n');
    ej2();
    fprintf('\n\n');
    fprintf('\n\n');
    ej3();
    fprintf('\n\n');
    fprintf('\n\n');
    ej4();
    fprintf('\n\n');
    fprintf('\n\n');
    ej5();
    fprintf('\n\n');
    fprintf('\n\n');
    ej6();
    fprintf('\n\n');
    fprintf('\n\n');
    ej7();
end

% Función para implementar la Regla del Trapecio
function integral_aproximada = regla_trapecio(f, a, b, n)
    % f: función a integrar
    % a: límite inferior
    % b: límite superior
    % n: número de subdivisiones
    h = (b-a)/n;
    suma = f(a) + f(b);
    for i=1:n-1
      suma = suma+2*f(a+i*h);
    end
      integral_aproximada = (h/2) * suma;
end

% Ejercicio 2 Método de Simpson Múltiple
function ej2()
    format long;  % 8 decimales

    % Función a integrar: sin(πt)
    a = 0;
    b = 0.5;
    integral_exacta = 1/pi;
    valores_n = [4, 6, 8, 12];

    fprintf('=== EJERCICIO 2: REGLA DE SIMPSON MÚLTIPLE ===\n\n');
    fprintf('Función: sin(πt)\n');
    fprintf('Intervalo: [%.1f, %.1f]\n', a, b);
    fprintf('Integral exacta: 1/π = %.8f\n\n', integral_exacta);

    % Tabla de resultados
    fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', 'N', 'Δt', 'ISim', 'Error', 'Log(Δt)', 'Log(Er)', 'Pendiente');
    fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', '---', '-----', '-------', '-----', '-------', '-------', '--------');

    % Vectores para análisis de convergencia
    delta_t = zeros(1, length(valores_n));
    errores = zeros(1, length(valores_n));

    for i = 1:length(valores_n)
        N = valores_n(i);
        dt = (b-a)/N;

        % Discretización del intervalo
        tg = linspace(a, b, N+1);
        yg = sin(pi*tg);

        % Algoritmo de Simpson Múltiple
        ISim = 0;
        for k = 2:2:N
            ISim = ISim + dt*(yg(k-1) + 4*yg(k) + yg(k+1))/3;
        end

        % Calcular error
        error_simpson = abs(ISim - integral_exacta);

        % Almacenar para análisis de convergencia
        delta_t(i) = dt;
        errores(i) = error_simpson;

        % Mostrar resultados en tabla
        log_dt = log10(dt);
        log_er = log10(error_simpson);

        if i > 1
            pendiente = (log_er - log10(errores(i-1))) / (log_dt - log10(delta_t(i-1)));
            fprintf('%-5d %-8.5f %-15.8f %-12.6e %-10.5f %-10.5f %-10.6f\n', N, dt, ISim, error_simpson, log_dt, log_er, pendiente);
        else
            fprintf('%-5d %-8.5f %-15.8f %-12.6e %-10.5f %-10.5f %-10s\n', N, dt, ISim, error_simpson, log_dt, log_er, '---');
        end
    end
end

% Función para implementar la Suma de Riemann
function integral_aproximada = suma_riemann(f, a, b, n)
    % f: función a integrar
    % a: límite inferior
    % b: límite superior
    % n: número de subdivisiones
    h = (b-a)/n;
    suma = 0;
    for i = 0:n-1  % Suma de Riemann usa N puntos (0 a N-1)
        suma = suma + f(a + i*h);
    end
    integral_aproximada = h * suma;
end

% Ejercicio 3 - Integración numérica de función exponencial
function ej3()
    fprintf('=== EJERCICIO 3: INTEGRACIÓN DE FUNCIÓN EXPONENCIAL ===\n\n');

    % Integral: ∫[0,6/w] e^(-wt) dt = (1/w)(1-e^(-6))

    valores_w = [2, 0.5];  % Dos casos de estudio
    valores_n = [6, 12, 18, 24];

    % Matrices para almacenar resultados de comparación
    resultados_riemann = zeros(length(valores_w), length(valores_n));
    resultados_trapecio = zeros(length(valores_w), length(valores_n));
    errores_riemann = zeros(length(valores_w), length(valores_n));
    errores_trapecio = zeros(length(valores_w), length(valores_n));

    % Almacenar valores exactos para cada w
    valores_exactos = zeros(1, length(valores_w));
    for j = 1:length(valores_w)
      w = valores_w(j);
      fprintf('=== CASO w = %.1f ===\n', w);
      f = @(t) exp(-w*t);
      a = 0;
      b = 6/w;
      integral_exacta = (1/w)*(1-exp(-6));
      valores_exactos(j) = integral_exacta;

      fprintf('Función: e^(-%.1ft), Intervalo: [%.1f, %.1f], Integral exacta: %.6f\n\n', w, a, b, integral_exacta);

      % TABLA 1: SUMA DE RIEMANN
      fprintf('TABLA 1: SUMA DE RIEMANN (w=%.1f)\n', w);
      fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', 'N', 'Δt', 'Im(Δt)', 'Error', 'Log(Δt)', 'Log(Er)', 'Pendiente');
      fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', '---', '-----', '-------', '-----', '-------', '-------', '--------');

      % Vectores temporales para cada w
      delta_t_temp = zeros(1, length(valores_n));
      errores_temp = zeros(1, length(valores_n));

      for i = 1:length(valores_n)
          n = valores_n(i);
          dt = (b-a)/n;  % Δt
          resultado_riemann = suma_riemann(f, a, b, n);
          error_riemann = abs(resultado_riemann - integral_exacta);

          % Almacenar en matrices de comparación
          resultados_riemann(j, i) = resultado_riemann;
          errores_riemann(j, i) = error_riemann;

          % Almacenar en vectores temporales para análisis
          delta_t_temp(i) = dt;
          errores_temp(i) = error_riemann;

          % Calcular logaritmos
          log_dt = log10(dt);
          log_er = log10(error_riemann);

          % Calcular pendiente (orden de error) con el punto anterior
          if i > 1
              pendiente = (log_er - log10(errores_temp(i-1))) / (log_dt - log10(delta_t_temp(i-1)));
              fprintf('%-5d %-8.5f %-15.8f %-12.6f %-10.5f %-10.5f %-10.6f\n', n, dt, resultado_riemann, error_riemann, log_dt, log_er, pendiente);
          else
              fprintf('%-5d %-8.5f %-15.8f %-12.6f %-10.5f %-10.5f %-10s\n', n, dt, resultado_riemann, error_riemann, log_dt, log_er, '---');
          end
      end

      fprintf('\n');

      % TABLA 2: TRAPECIOS MÚLTIPLE (ITrap)
      fprintf('TABLA 2: TRAPECIOS MÚLTIPLE (w=%.1f)\n', w);
      fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', 'N', 'Δt', 'ITrap(Δt)', 'Error', 'Log(Δt)', 'Log(Er)', 'Pendiente');
      fprintf('%-5s %-8s %-15s %-12s %-10s %-10s %-10s\n', '---', '-----', '-------', '-----', '-------', '-------', '--------');

      % Reinicializar vectores temporales para Trapecios
      delta_t_temp = zeros(1, length(valores_n));
      errores_temp = zeros(1, length(valores_n));

      for i = 1:length(valores_n)
          n = valores_n(i);
          dt = (b-a)/n;  % Δt
          resultado_trapecio = regla_trapecio(f, a, b, n);
          error_trapecio = abs(resultado_trapecio - integral_exacta);

          % Almacenar en matrices de comparación
          resultados_trapecio(j, i) = resultado_trapecio;
          errores_trapecio(j, i) = error_trapecio;

          % Almacenar en vectores temporales para análisis
          delta_t_temp(i) = dt;
          errores_temp(i) = error_trapecio;

          % Calcular logaritmos
          log_dt = log10(dt);
          log_er = log10(error_trapecio);

          % Calcular pendiente (orden de error) con el punto anterior
          if i > 1
              pendiente = (log_er - log10(errores_temp(i-1))) / (log_dt - log10(delta_t_temp(i-1)));
              fprintf('%-5d %-8.5f %-15.8f %-12.6f %-10.5f %-10.5f %-10.6f\n', n, dt, resultado_trapecio, error_trapecio, log_dt, log_er, pendiente);
          else
              fprintf('%-5d %-8.5f %-15.8f %-12.6f %-10.5f %-10.5f %-10s\n', n, dt, resultado_trapecio, error_trapecio, log_dt, log_er, '---');
          end
      end

      fprintf('\n\n');
    end

    % ========== COMPARACIÓN Y GRÁFICOS ==========
    fprintf('=== ANÁLISIS COMPARATIVO ===\n');
    crear_graficos_comparacion(valores_w, valores_n, resultados_riemann, resultados_trapecio, valores_exactos);
end

% Función para crear gráficos comparativos del Ejercicio 3
function crear_graficos_comparacion(valores_w, valores_n, resultados_riemann, resultados_trapecio, valores_exactos)
    % Configurar toolkit gráfico recomendado por Octave
    % graphics_toolkit('qt');  % Comentado - usar toolkit por defecto disponible

    % Crear figura con subplots para comparación
    figure('Name', 'Ejercicio 3: Comparación de Métodos vs w', 'Position', [100, 100, 1200, 800]);

    % Subplot 1: Convergencia para w=2
    subplot(2,2,1);
    plot(valores_n, resultados_trapecio(1,:), 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    plot(valores_n, resultados_riemann(1,:), 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
    plot([min(valores_n), max(valores_n)], [valores_exactos(1), valores_exactos(1)], 'k--', 'LineWidth', 1.5);
    title('Convergencia para w = 2.0');
    xlabel('Número de intervalos (N)');
    ylabel('Valor de la integral');
    legend('Trapecios Múltiple', 'Suma de Riemann', 'Location', 'northeast');
    grid on;

    % Subplot 2: Convergencia para w=0.5
    subplot(2,2,2);
    plot(valores_n, resultados_trapecio(2,:), 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    plot(valores_n, resultados_riemann(2,:), 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
    plot([min(valores_n), max(valores_n)], [valores_exactos(2), valores_exactos(2)], 'k--', 'LineWidth', 1.5);
    title('Convergencia para w = 0.5');
    xlabel('Número de intervalos (N)');
    ylabel('Valor de la integral');
    legend('Trapecios Múltiple', 'Suma de Riemann', 'Location', 'northeast');
    grid on;

    % Subplot 3: Comparación de errores para Trapecios
    subplot(2,2,3);
    errores_trap_w2 = abs(resultados_trapecio(1,:) - valores_exactos(1));
    errores_trap_w05 = abs(resultados_trapecio(2,:) - valores_exactos(2));
    semilogy(valores_n, errores_trap_w2, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    semilogy(valores_n, errores_trap_w05, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    title('Errores - Método de Trapecios');
    xlabel('Número de intervalos (N)');
    ylabel('Error absoluto (escala log)');
    legend('w = 2.0', 'w = 0.5', 'Location', 'northeast');
    grid on;

    % Subplot 4: Comparación de errores para Riemann
    subplot(2,2,4);
    errores_riem_w2 = abs(resultados_riemann(1,:) - valores_exactos(1));
    errores_riem_w05 = abs(resultados_riemann(2,:) - valores_exactos(2));
    semilogy(valores_n, errores_riem_w2, 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    semilogy(valores_n, errores_riem_w05, 'm-d', 'LineWidth', 2, 'MarkerSize', 8);
    title('Errores - Suma de Riemann');
    xlabel('Número de intervalos (N)');
    ylabel('Error absoluto (escala log)');
    legend('w = 2.0', 'w = 0.5', 'Location', 'northeast');
    grid on;

    % Análisis textual de la comparación
    fprintf('CONCLUSIONES DEL ANÁLISIS COMPARATIVO:\n');
    fprintf('=====================================\n\n');

    % Comparar convergencia entre métodos
    error_final_trap_w2 = errores_trap_w2(end);
    error_final_riem_w2 = errores_riem_w2(end);
    error_final_trap_w05 = errores_trap_w05(end);
    error_final_riem_w05 = errores_riem_w05(end);

    fprintf('1. PRECISIÓN FINAL (N=24):\n');
    fprintf('   w=2.0:  Trapecios = %.2e, Riemann = %.2e\n', error_final_trap_w2, error_final_riem_w2);
    fprintf('   w=0.5:  Trapecios = %.2e, Riemann = %.2e\n\n', error_final_trap_w05, error_final_riem_w05);

    fprintf('2. IMPACTO DEL PARÁMETRO w:\n');
    if error_final_trap_w2 < error_final_trap_w05
        fprintf('   - Trapecios converge MEJOR con w=2.0 (decaimiento rápido)\n');
    else
        fprintf('   - Trapecios converge MEJOR con w=0.5 (decaimiento lento)\n');
    end

    if error_final_riem_w2 < error_final_riem_w05
        fprintf('   - Riemann converge MEJOR con w=2.0 (decaimiento rápido)\n\n');
    else
        fprintf('   - Riemann converge MEJOR con w=0.5 (decaimiento lento)\n\n');
    end

    fprintf('3. COMPARACIÓN ENTRE MÉTODOS:\n');
    if error_final_trap_w2 < error_final_riem_w2
        fprintf('   - Para w=2.0: TRAPECIOS es más preciso\n');
    else
        fprintf('   - Para w=2.0: RIEMANN es más preciso\n');
    end

    if error_final_trap_w05 < error_final_riem_w05
        fprintf('   - Para w=0.5: TRAPECIOS es más preciso\n\n');
    else
        fprintf('   - Para w=0.5: RIEMANN es más preciso\n\n');
    end
end

% Ejercicio 4 - Integración de funciones trigonométricas - Ortogonalidad
function ej4()
    fprintf('=== EJERCICIO 4: INTEGRACIÓN DE FUNCIONES TRIGONOMÉTRICAS ===\n\n');
    fprintf('Verificación de propiedades de ortogonalidad\n\n');

    % Valores de N según el PDF (parte c): "N=10; 20; 30; 40"
    valores_n = [10, 20, 30, 40];

    % Casos de prueba específicos del PDF
    casos_ortogonalidad = [
        % [k, m, tipo] donde tipo: 1=sin*cos, 2=sin*sin, 3=cos*cos
        1, 1, 1;   % sin(θ)*cos(θ) - debe dar 0
        1, 2, 1;   % sin(θ)*cos(2θ) - debe dar 0
        1, 1, 2;   % sin(θ)*sin(θ) - debe dar π
        2, 2, 2;   % sin(2θ)*sin(2θ) - debe dar π
        1, 2, 2;   % sin(θ)*sin(2θ) - debe dar 0
        0, 0, 3;   % cos(0)*cos(0) = 1*1 - debe dar 2π
        1, 1, 3;   % cos(θ)*cos(θ) - debe dar π
        1, 2, 3;   % cos(θ)*cos(2θ) - debe dar 0
    ];

    fprintf('Casos de verificación de ortogonalidad:\n');
    fprintf('=====================================\n\n');

    % TODO(human): Implementa la verificación de ortogonalidad según el PDF
    % Para cada caso en casos_ortogonalidad:
    %   1. Extrae k, m, tipo: k=casos_ortogonalidad(i,1), m=casos_ortogonalidad(i,2), tipo=casos_ortogonalidad(i,3)
    %   2. Define función según tipo: @(theta) sin(k*theta).*cos(m*theta), etc.
    %   3. Calcula valor teórico según reglas del PDF
    %   4. Para cada N: aplica regla_trapecio y suma_riemann en [0, 2π]
    %   5. Muestra tabla comparativa con errores y analiza convergencia
    for i = 1: length(casos_ortogonalidad)
      k = casos_ortogonalidad(i, 1);
      m = casos_ortogonalidad(i, 2);
      t = casos_ortogonalidad(i, 3);
      % Crea la función según el tipo
      if t == 1
        f = @(theta) sin(k*theta) .* cos(m*theta);
      elseif t == 2
        f = @(theta) sin(k*theta) .* sin(m*theta);
      elseif t == 3
        f = @(theta) cos(k*theta) .* cos(m*theta);
      end

      % Calcular valor teórico
      if t == 1
        valor_teorico = 0;
      elseif t == 2
        if k == 0 && m == 0
          valor_teorico = 0;
        elseif k == m && k ~= 0
          valor_teorico = (2*pi)/2;
        else
          valor_teorico = 0;
        end
      elseif t == 3
        if k == 0 && m == 0
          valor_teorico = 2*pi;
        elseif k == m && k ~= 0
          valor_teorico = (2*pi)/2;
        else
          valor_teorico = 0;
        end
      end
      % Mostrar información del caso
      fprintf('Caso %d: ', i);
      if t == 1
          fprintf('sin(%dθ) × cos(%dθ) = %.4f (teórico)\n', k, m, valor_teorico);
      elseif t == 2
          fprintf('sin(%dθ) × sin(%dθ) = %.4f (teórico)\n', k, m, valor_teorico);
      elseif t == 3
          fprintf('cos(%dθ) × cos(%dθ) = %.4f (teórico)\n', k, m, valor_teorico);
      end
      for j = 1:length(valores_n)
        N = valores_n(j);

        % Aplicando métodos en intervalos [0, 2π]
        resultado_trapecio = regla_trapecio(f, 0, 2*pi, N);
        resultado_riemann = suma_riemann(f, 0, 2*pi, N);

        % Calcular errores
        error_trapecio = abs(resultado_trapecio - valor_teorico);
        error_riemann = abs(resultado_riemann - valor_teorico);

        % Mostrar resultados
        fprintf('  N=%2d: Trapecio=%.6f (err=%.2e), Riemann=%.6f (err=%.2e)\n', ...
                N, resultado_trapecio, error_trapecio, resultado_riemann, error_riemann);
      end
      fprintf('\n');
    end
end

% Ejercicio 5 - Obtención de la función integral de sen(pi*t)
function ej5()
    format long;  % 8 decimales según especificación

    fprintf('=== EJERCICIO 5: OBTENCIÓN DE LA FUNCIÓN INTEGRAL DE sen(π t) ===\n\n');
    fprintf('Búsqueda de la función primitiva G(t) tal que dG(t)/dt = sen(π t)\n');
    fprintf('Condición inicial: G(0) = -1/π\n\n');

    %% PARTE A: Verificación de equivalencia con N=5 usando datos exactos del PDF
    fprintf('=== PARTE A: Verificación de equivalencia con N=5 ===\n');
    N = 5;
    dt = 0.1;  % Δt = 0.5/5 = 0.1

    % DATOS EXACTOS del PDF (página 7)
    tg = [0.00000, 0.10000, 0.20000, 0.30000, 0.40000, 0.50000];
    yg = [0.00000, 0.30902, 0.58779, 0.80902, 0.95106, 1.00000];

    fprintf('Datos exactos del PDF:\n');
    fprintf('tg = [%.5f %.5f %.5f %.5f %.5f %.5f]\n', tg);
    fprintf('yg = [%.5f %.5f %.5f %.5f %.5f %.5f]\n\n', yg);

    % Método 1: Formulación matricial EXACTA del PDF
    fprintf('--- Método 1: Formulación Matricial del PDF ---\n');
    Gd_matricial = calcular_primitiva_matricial(yg, dt, N);

    fprintf('Resultados método matricial:\n');
    for i = 1:length(Gd_matricial)
        fprintf('Gd(%d) = %.8f\n', i, Gd_matricial(i));
    end
    fprintf('\n');

    % Método 2: Algoritmo iterativo EXACTO del PDF
    fprintf('--- Método 2: Algoritmo Iterativo del PDF ---\n');
    Gd_iterativo = calcular_primitiva_iterativo(yg, dt, N);

    fprintf('Resultados método iterativo:\n');
    for i = 1:length(Gd_iterativo)
        fprintf('Gd(%d) = %.8f\n', i, Gd_iterativo(i));
    end
    fprintf('\n');

    % Verificar equivalencia (punto a del PDF)
    diferencia_maxima = max(abs(Gd_matricial - Gd_iterativo));
    fprintf('--- Verificación de Equivalencia (Punto a del PDF) ---\n');
    fprintf('Diferencia máxima entre métodos: %.2e\n', diferencia_maxima);
    if diferencia_maxima < 1e-12
        fprintf('✓ EQUIVALENCIA VERIFICADA: Los métodos son equivalentes\n\n');
    else
        fprintf('✗ ERROR: Los métodos NO son equivalentes\n');
        fprintf('Diferencias por elemento:\n');
        for i = 1:length(Gd_matricial)
            diff = abs(Gd_matricial(i) - Gd_iterativo(i));
            fprintf('  Gd(%d): Matricial = %.8f, Iterativo = %.8f, Diff = %.2e\n', ...
                    i, Gd_matricial(i), Gd_iterativo(i), diff);
        end
        fprintf('\n');
    end

    %% PARTE B: Visualización con N=200
    fprintf('=== PARTE B: Visualización Gráfica con N=200 ===\n');
    N_grafico = 200;
    dt_periodo = 1/N_grafico;  % Para período completo [0,1]

    % Para visualización completa de un período [0,1] como especifica el PDF
    t_periodo = linspace(0, 1, N_grafico+1);  % Período completo [0,1]
    y_periodo = sin(pi*t_periodo);

    % Calcular primitiva para período completo usando algoritmo iterativo del PDF
    G_periodo = calcular_primitiva_iterativo_PDF(y_periodo, dt_periodo, N_grafico);

    % Crear visualización gráfica
    crear_grafica_primitiva(t_periodo, y_periodo, G_periodo);

    fprintf('Gráfica generada mostrando:\n');
    fprintf('- Azul: función original sen(π t) en un período completo [0,1]\n');
    fprintf('- Rojo: función primitiva G(t)\n\n');
end
% Función para calcular primitiva usando formulación matricial EXACTA del PDF
function Gd = calcular_primitiva_matricial(yg, dt, N)
    % Implementa la formulación matricial EXACTA mostrada en la página 7 del PDF:
    % Gd = (-1/π)*[1;1;1;1;1;1] + (Δt/2) * matriz_triangular * yg

    % Vector condición inicial
    G0 = -1/pi;
    vector_unos = ones(N+1, 1);

    % Matriz triangular EXACTA del PDF (página 7):
    % ⎡0 0 0 0 0 0⎤
    % ⎢1 1 0 0 0 0⎥
    % ⎢1 2 1 0 0 0⎥
    % ⎢1 2 2 1 0 0⎥
    % ⎢1 2 2 2 1 0⎥
    % ⎣1 2 2 2 2 1⎦

    matriz_triangular = zeros(N+1, N+1);

    % Construir matriz triangular según el patrón del PDF
    for i = 1:N+1  % filas
        for j = 1:N+1  % columnas
            if i == 1
                matriz_triangular(i, j) = 0;  % Primera fila toda ceros
            elseif j == 1
                matriz_triangular(i, j) = 1;  % Primera columna (excepto fila 1) todos unos
            elseif j == i
                matriz_triangular(i, j) = 1;  % Diagonal todos unos
            elseif j < i
                matriz_triangular(i, j) = 2;  % Parte triangular inferior todos doses
            else
                matriz_triangular(i, j) = 0;  % Parte triangular superior todos ceros
            end
        end
    end

    % Aplicar la fórmula del PDF: Gd = G0*[1;1;...] + (dt/2)*matriz*yg
    yg_column = yg(:);  % Asegurar que yg es vector columna
    Gd = G0 * vector_unos + (dt/2) * matriz_triangular * yg_column;
end

% Función para calcular primitiva usando algoritmo iterativo EXACTO del PDF
function Gd = calcular_primitiva_iterativo(yg, dt, N)
    % Algoritmo EXACTO según especificación del PDF (página 7):
    % Gd(1) = -1/π; Trap=0;
    % for k=1:N
    %   Trap = dt*(yg(k)+yg(k+1))/2;
    %   Gd(k+1) = Gd(k) + Trap
    % end
    % Nota: El PDF tiene un error tipográfico (dice -1/w pero debe ser -1/π)

    Gd = zeros(N+1, 1);
    Gd(1) = -1/pi;  % Condición inicial (corrigiendo error tipográfico del PDF)
    Trap = 0;

    for k = 1:N
        Trap = dt * (yg(k) + yg(k+1)) / 2;
        Gd(k+1) = Gd(k) + Trap;
    end
end

%Función para calcular primitiva usando formulación matricial EXACTA del PDF
function Gd = calcular_primitiva_matricial_PDF(yg, dt, N)
    % Implementa la formulación matricial EXACTA mostrada en la página 7 del PDF:
    % Gd = (-1/π)*[1;1;1;1;1;1] + (Δt/2) * matriz_triangular * yg

    % Vector condición inicial
    G0 = -1/pi;
    vector_unos = ones(N+1, 1);

    % Matriz triangular EXACTA del PDF (página 7):
    % ⎡0 0 0 0 0 0⎤
    % ⎢1 1 0 0 0 0⎥
    % ⎢1 2 1 0 0 0⎥
    % ⎢1 2 2 1 0 0⎥
    % ⎢1 2 2 2 1 0⎥
    % ⎣1 2 2 2 2 1⎦

    matriz_triangular = zeros(N+1, N+1);

    % Construir matriz triangular según el patrón del PDF
    for i = 1:N+1  % filas
        for j = 1:N+1  % columnas
            if i == 1
                matriz_triangular(i, j) = 0;  % Primera fila toda ceros
            elseif j == 1
                matriz_triangular(i, j) = 1;  % Primera columna (excepto fila 1) todos unos
            elseif j == i
                matriz_triangular(i, j) = 1;  % Diagonal todos unos
            elseif j < i
                matriz_triangular(i, j) = 2;  % Parte triangular inferior todos doses
            else
                matriz_triangular(i, j) = 0;  % Parte triangular superior todos ceros
            end
        end
    end

    % Aplicar la fórmula del PDF: Gd = G0*[1;1;...] + (dt/2)*matriz*yg
    yg_column = yg(:);  % Asegurar que yg es vector columna
    Gd = G0 * vector_unos + (dt/2) * matriz_triangular * yg_column;
end

% Función para calcular primitiva usando algoritmo iterativo EXACTO del PDF
function Gd = calcular_primitiva_iterativo_PDF(yg, dt, N)
    % Algoritmo EXACTO según especificación del PDF (página 7):
    % Gd(1) = -1/π; Trap=0;
    % for k=1:N
    %   Trap = dt*(yg(k)+yg(k+1))/2;
    %   Gd(k+1) = Gd(k) + Trap
    % end
    % Nota: El PDF tiene un error tipográfico (dice -1/w pero debe ser -1/π)

    Gd = zeros(N+1, 1);
    Gd(1) = -1/pi;  % Condición inicial (corrigiendo error tipográfico del PDF)
    Trap = 0;

    for k = 1:N
        Trap = dt * (yg(k) + yg(k+1)) / 2;
        Gd(k+1) = Gd(k) + Trap;
    end
end


% Función para crear la gráfica comparativa
function crear_grafica_primitiva(t, y_funcion, G_primitiva)
    figure('Name', 'Ejercicio 5: Función sen(π t) y su Primitiva', 'Position', [100, 100, 1000, 600]);

    plot(t, y_funcion, 'b-', 'LineWidth', 2, 'DisplayName', 'sen(π t)');
    hold on;
    plot(t, G_primitiva, 'r-', 'LineWidth', 2, 'DisplayName', 'G(t) - Primitiva');

    title('Función sen(π t) y su Función Primitiva G(t)');
    xlabel('tiempo t');
    ylabel('Amplitud');
    legend('Location', 'northeast');
    grid on;

    % Configurar ejes para mejor visualización
    axis([0 1 -1.5 1.5]);

    fprintf('Gráfica generada mostrando:\n');
    fprintf('- Azul: función original sen(π t)\n');
    fprintf('- Rojo: función primitiva G(t)\n\n');
end

% Ejercicio 6 - Obtención de la función integral de y(t) = e^(-pt)
function ej6()
    format long;

    fprintf('=== EJERCICIO 6: OBTENCIÓN DE LA FUNCIÓN INTEGRAL DE y(t) = e^(-pt) ===\n\n');
    fprintf('Función: y(t) = e^(-pt)\n');
    fprintf('Primitiva analítica: G(t) = (1/p)(1 - e^(-pt)), G(0) = 0\n\n');

    %% CASO 1: p = 2, N = 10
    fprintf('=== CASO 1: p = 2, N = 10 ===\n');
    p1 = 2;
    N1 = 10;
    t_final1 = 3;  % Dominio [0,3] para capturar ~99% del decaimiento
    dt1 = t_final1 / N1;

    % Discretización
    t1 = linspace(0, t_final1, N1+1);
    y1 = exp(-p1 * t1);

    % Calcular primitiva numérica
    G1_numerica = calcular_primitiva_iterativo_PDF(y1, dt1, N1);
    G1_numerica(1) = 0;  % Condición inicial G(0) = 0 para este ejercicio

    % Primitiva analítica
    G1_analitica = (1/p1) * (1 - exp(-p1 * t1));

    % Generar gráfica
    crear_grafica_ejercicio6(t1, y1, G1_numerica, G1_analitica, p1, N1, 1);

    %% CASO 2: p = 0.5, N = 10
    fprintf('=== CASO 2: p = 0.5, N = 10 ===\n');
    p2 = 0.5;
    N2 = 10;
    t_final2 = 12;  % Dominio [0,12] para decaimiento lento
    dt2 = t_final2 / N2;

    % Discretización
    t2 = linspace(0, t_final2, N2+1);
    y2 = exp(-p2 * t2);

    % Calcular primitiva numérica
    G2_numerica = calcular_primitiva_iterativo_PDF(y2, dt2, N2);
    G2_numerica(1) = 0;  % Condición inicial G(0) = 0

    % Primitiva analítica
    G2_analitica = (1/p2) * (1 - exp(-p2 * t2));

    % Valor teórico en estado estacionario
    valor_teorico = 1/p2;  % = 2
    valor_numerico_final = G2_numerica(end);

    fprintf('Valor teórico en estado estacionario: %.6f\n', valor_teorico);
    fprintf('Valor numérico final (N=10): %.6f\n', valor_numerico_final);
    fprintf('Error absoluto: %.6f\n', abs(valor_numerico_final - valor_teorico));
    fprintf('Como se esperaba, el valor numérico es SUPERIOR al teórico.\n\n');

    % Generar gráfica
    crear_grafica_ejercicio6(t2, y2, G2_numerica, G2_analitica, p2, N2, 2);

    %% CASO 3: p = 0.5, N = 20
    fprintf('=== CASO 3: p = 0.5, N = 20 ===\n');
    p3 = 0.5;
    N3 = 20;
    t_final3 = 12;  % Mismo dominio
    dt3 = t_final3 / N3;

    % Discretización
    t3 = linspace(0, t_final3, N3+1);
    y3 = exp(-p3 * t3);

    % Calcular primitiva numérica
    G3_numerica = calcular_primitiva_iterativo_PDF(y3, dt3, N3);
    G3_numerica(1) = 0;  % Condición inicial G(0) = 0

    % Primitiva analítica
    G3_analitica = (1/p3) * (1 - exp(-p3 * t3));

    % Valor teórico en estado estacionario
    valor_numerico_final_N20 = G3_numerica(end);

    fprintf('Valor teórico en estado estacionario: %.6f\n', valor_teorico);
    fprintf('Valor numérico final (N=20): %.6f\n', valor_numerico_final_N20);
    fprintf('Error absoluto: %.6f\n', abs(valor_numerico_final_N20 - valor_teorico));
    fprintf('Como se esperaba, la diferencia con el valor teórico DISMINUYE.\n\n');

    % Generar gráfica
    crear_grafica_ejercicio6(t3, y3, G3_numerica, G3_analitica, p3, N3, 3);

    %% COMPARACIÓN DE CONVERGENCIA
    fprintf('=== COMPARACIÓN DE CONVERGENCIA ===\n');
    fprintf('Error N=10:  %.6f\n', abs(valor_numerico_final - valor_teorico));
    fprintf('Error N=20:  %.6f\n', abs(valor_numerico_final_N20 - valor_teorico));
    fprintf('Mejora:     %.2fx\n', abs(valor_numerico_final - valor_teorico) / abs(valor_numerico_final_N20 - valor_teorico));
end

% Función auxiliar para crear gráficas del ejercicio 6
function crear_grafica_ejercicio6(t, y_funcion, G_numerica, G_analitica, p, N, caso_num)
    figure('Name', sprintf('Ejercicio 6 - Caso %d: p=%.1f, N=%d', caso_num, p, N), ...
           'Position', [100 + caso_num*50, 100 + caso_num*50, 1000, 600]);

    subplot(2,1,1);
    plot(t, y_funcion, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('y(t) = e^{-%.1ft}', p));
    title(sprintf('Función Original: y(t) = e^{-%.1ft}', p));
    xlabel('tiempo t');
    ylabel('y(t)');
    legend('Location', 'northeast');
    grid on;

    subplot(2,1,2);
    plot(t, G_numerica, 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', sprintf('G(t) Numérica N=%d', N));
    hold on;
    plot(t, G_analitica, 'g--', 'LineWidth', 2, 'DisplayName', 'G(t) Analítica');

    % Marcar valor en estado estacionario
    valor_teorico = 1/p;
    plot([0, t(end)], [valor_teorico, valor_teorico], 'k:', 'LineWidth', 1, ...
         'DisplayName', sprintf('Valor teórico = %.1f', valor_teorico));

    title(sprintf('Primitiva G(t): p=%.1f, N=%d', p, N));
    xlabel('tiempo t');
    ylabel('G(t)');
    legend('Location', 'southeast');
    grid on;

    fprintf('Gráfica generada para caso %d: p=%.1f, N=%d\n', caso_num, p, N);
end

% Ejercicio 7 - Obtención de la función integral de y(t) = e^(-pt)sin(wt)
function ej7()
    format long;

    fprintf('=== EJERCICIO 7: OBTENCIÓN DE LA FUNCIÓN INTEGRAL DE y(t) = e^(-pt)sin(wt) ===\n\n');
    fprintf('Función: y(t) = e^(-pt)sin(wt)\n');
    fprintf('Parámetros: p = 0.5, w = 4\n');
    fprintf('Dominio: t ∈ [0, 8]\n');
    fprintf('Condición inicial: G(0) = 0\n\n');

    % Parámetros del problema
    p = 0.5;
    w = 4;
    t_final = 8;

    %% CASO 1: N = 50
    fprintf('=== CASO 1: N = 50 ===\n');
    N1 = 50;
    dt1 = t_final / N1;

    % Discretización
    t1 = linspace(0, t_final, N1+1);
    y1 = exp(-p * t1) .* sin(w * t1);

    % Calcular primitiva numérica
    G1_numerica = calcular_primitiva_iterativo_PDF(y1, dt1, N1);
    G1_numerica(1) = 0;  % Condición inicial G(0) = 0

    fprintf('Primitiva calculada con N=50\n');
    fprintf('G(0) = %.8f\n', G1_numerica(1));
    fprintf('G(8) = %.8f\n', G1_numerica(end));
    fprintf('\n');

    %% CASO 2: N = 100
    fprintf('=== CASO 2: N = 100 ===\n');
    N2 = 100;
    dt2 = t_final / N2;

    % Discretización
    t2 = linspace(0, t_final, N2+1);
    y2 = exp(-p * t2) .* sin(w * t2);

    % Calcular primitiva numérica
    G2_numerica = calcular_primitiva_iterativo_PDF(y2, dt2, N2);
    G2_numerica(1) = 0;  % Condición inicial G(0) = 0

    fprintf('Primitiva calculada con N=100\n');
    fprintf('G(0) = %.8f\n', G2_numerica(1));
    fprintf('G(8) = %.8f\n', G2_numerica(end));
    fprintf('\n');

    %% ANÁLISIS DE DIFERENCIAS CERCA DE t=0
    fprintf('=== ANÁLISIS DE DIFERENCIAS CERCA DE t=0 ===\n');

    % Interpolar G1 en los puntos de t2 para comparar
    G1_interpolada = interp1(t1, G1_numerica, t2, 'linear');

    % Calcular diferencias
    diferencias = abs(G2_numerica - G1_interpolada);

    % Encontrar diferencia máxima cerca de t=0 (primeros 10% del dominio)
    indices_cercanos_a_0 = t2 <= 0.1 * t_final;
    diferencia_maxima_t0 = max(diferencias(indices_cercanos_a_0));
    [~, idx_max] = max(diferencias(indices_cercanos_a_0));
    t_max_diff = t2(indices_cercanos_a_0);
    t_max_diff = t_max_diff(idx_max);

    fprintf('Diferencia máxima cerca de t=0: %.8f en t=%.4f\n', diferencia_maxima_t0, t_max_diff);
    fprintf('Diferencia máxima total: %.8f\n', max(diferencias));

    % Mostrar primeros valores para verificar diferencias cerca de t=0
    fprintf('\nPrimeros valores (cerca de t=0):\n');
    fprintf('t      | G_N50     | G_N100    | Diferencia\n');
    fprintf('-------|-----------|-----------|----------\n');
    for i = 1:min(6, length(t2))
        fprintf('%.3f  | %.6f  | %.6f  | %.6f\n', ...
                t2(i), G1_interpolada(i), G2_numerica(i), diferencias(i));
    end
    fprintf('\n');

    %% GENERAR GRÁFICAS COMPARATIVAS
    crear_grafica_ejercicio7(t1, y1, G1_numerica, t2, y2, G2_numerica, p, w, N1, N2);

    fprintf('Como se observa en las gráficas, existen diferencias notables\n');
    fprintf('entre ambas aproximaciones, especialmente para valores de t cercanos a 0.\n');
    fprintf('Esto se debe a que la alta frecuencia (w=4) requiere mayor resolución temporal.\n\n');
end

% Función auxiliar para crear gráficas comparativas del ejercicio 7
function crear_grafica_ejercicio7(t1, y1, G1, t2, y2, G2, p, w, N1, N2)
    % Crear figura con subplots side-by-side como en el PDF
    figure('Name', sprintf('Ejercicio 7: e^{-%.1ft}sin(%.0ft) - Comparación N=%d vs N=%d', p, w, N1, N2), ...
           'Position', [100, 100, 1400, 800]);

    % Subplot para N=50
    subplot(2,2,1);
    plot(t1, y1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('y(t) = e^{-%.1ft}sin(%.0ft)', p, w));
    title(sprintf('Función Original - N=%d', N1));
    xlabel('tiempo t');
    ylabel('y(t)');
    legend('Location', 'northeast');
    grid on;
    axis([0 8 -1 1]);

    subplot(2,2,3);
    plot(t1, y1, 'b-', 'LineWidth', 1, 'DisplayName', sprintf('y(t) = e^{-%.1ft}sin(%.0ft)', p, w));
    hold on;
    plot(t1, G1, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('G(t) Numérica N=%d', N1));
    title(sprintf('Primitiva G(t) - N=%d', N1));
    xlabel('tiempo t');
    ylabel('Amplitud');
    legend('Location', 'northeast');
    grid on;
    axis([0 8 -1 1]);

    % Subplot para N=100
    subplot(2,2,2);
    plot(t2, y2, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('y(t) = e^{-%.1ft}sin(%.0ft)', p, w));
    title(sprintf('Función Original - N=%d', N2));
    xlabel('tiempo t');
    ylabel('y(t)');
    legend('Location', 'northeast');
    grid on;
    axis([0 8 -1 1]);

    subplot(2,2,4);
    plot(t2, y2, 'b-', 'LineWidth', 1, 'DisplayName', sprintf('y(t) = e^{-%.1ft}sin(%.0ft)', p, w));
    hold on;
    plot(t2, G2, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('G(t) Numérica N=%d', N2));
    title(sprintf('Primitiva G(t) - N=%d', N2));
    xlabel('tiempo t');
    ylabel('Amplitud');
    legend('Location', 'northeast');
    grid on;
    axis([0 8 -1 1]);

    fprintf('Gráfica generada mostrando comparación N=%d vs N=%d\n', N1, N2);

    % Crear figura adicional para mostrar diferencias cerca de t=0
    figure('Name', 'Ejercicio 7: Diferencias cerca de t=0', ...
           'Position', [200, 200, 1000, 600]);

    % Interpolar para comparar
    G1_interpolada = interp1(t1, G1, t2, 'linear');
    diferencias = abs(G2 - G1_interpolada);

    subplot(2,1,1);
    plot(t2, G1_interpolada, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('G(t) N=%d (interpolada)', N1));
    hold on;
    plot(t2, G2, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('G(t) N=%d', N2));
    title('Comparación de Primitivas');
    xlabel('tiempo t');
    ylabel('G(t)');
    legend('Location', 'northeast');
    grid on;
    xlim([0 2]);  % Foco en los primeros valores

    subplot(2,1,2);
    % Asegurar que todas las diferencias sean positivas para el plot logarítmico
    diferencias_positivas = max(diferencias, 1e-16);  % Valor mínimo para evitar problemas con log
    semilogy(t2, diferencias_positivas, 'k-', 'LineWidth', 2, 'DisplayName', '|G_{N=100} - G_{N=50}|');
    title('Diferencia Absoluta entre Aproximaciones (escala logarítmica)');
    xlabel('tiempo t');
    ylabel('Diferencia absoluta');
    legend('Location', 'northeast');
    grid on;
    xlim([0 2]);  % Foco en los primeros valores donde las diferencias son más notables

    fprintf('Gráfica adicional generada mostrando diferencias cerca de t=0\n');
end

