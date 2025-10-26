function tp2 
    ej1();
    fprintf('\n\n');
    ej2();
    fprintf('\n\n');
    ej3();
    fprintf('\n\n');
    ej4();
    fprintf('\n\n');
endfunction

function ej1 
    % Función a analizar: y = sen(x)
    % Punto de evaluación x = pi /10
    % Valor primera derivada = cos(pi/10) = 0.951
    x0 = pi/10;
    f = @(x) sin(x);
    df_exact = @(x) cos(x);
    d2f_exact = @(x) -sin(x);

    %Valor de la primera derivada en x0 es:
    derivada_exacta_1 = df_exact(x0);
    fprintf('Valor de la primera derivada en x=π/10: %.9f\n', derivada_exacta_1);
    %Valor de la segunda derivada en x0 es:
    derivada_exacta_2 = d2f_exact(x0);
    fprintf('Valor de la segunda derivada en x=π/10: %.9f\n', derivada_exacta_2);

    % a) Derivda primera con dos puntos y su error asociado
    % FÓRMULAS: f'(x) = (f(x+h) - f(x)) / h
    % ERROR: E = -(h/2) - f''(e) = O(h)
    fprintf('\n ==== PARTE A: DERIVDA PRIMERA CON DOS VALORES ====\n');

    % Diferentes valores de h
    h_values=[0.00001, 0.0001, 0.001, 0.01, 0.1];

    %Inicializar valores para guardar los resultados
    deriv = zeros(size(h_values));
    error_derivada = zeros(size(h_values));

    fprintf('h\t\tDerivada Adelante\tError Adelante\t\tlog(h)\t\tlog(Error)\n');
    fprintf('---\t\t----------------\t---------------\t\t------\t\t----------\n');

    for i = 1 : length(h_values)
        h = h_values(i);

        % Calcular la primera derivda
        deriv(i) = (f(x0 + h) - f(x0)) / h;

        % Calcular error absoluto
        error_derivada(i) = abs(deriv(i) - derivada_exacta_1);

        % Mostrar resultados
        fprintf('%.5f\t%.9f\t\t%.2e\t\t%.1f\t\t%.3f\n', ...
            h, deriv(i), error_derivada(i), log10(h), log10(error_derivada(i)));
    end 

    % Calcular la pendiente de la recha log(error) vs log(h)
    log_h = log10(h_values);
    log_error = log10(error_derivada);
    pendiente_primera = polyfit(log_h, log_error, 1);
    fprintf('\nPendiente (orden de convergencia): %.2f\n', pendiente_primera(1));

    %% PARTE B) DIFERENCIAS CENTRALES (TRES PUNTOS) - PRIMERA DERIVADA
    % Fórmula: f'(x) ≈ (f(x+h) - f(x-h))/(2h)
    % Error: E = -(h²/6)*f'''(ξ) = O(h²)

    fprintf('\n=== PARTE B: DIFERENCIAS CENTRALES (PRIMERA DERIVADA) ===\n');
    derivada_central_1 = zeros(size(h_values));
    error_segunda_1 = zeros(size(h_values));

    fprintf('h\t\tDerivada Central\tError Central\t\tlog(h)\t\tlog(Error)\n');
    fprintf('---\t\t----------------\t---------------\t\t------\t\t----------\n');

    for i = 1:length(h_values)
        h = h_values(i);

        % Calcular segunda derivada central 
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
    h_values_2 = [0.1, 0.01, 0.001];

    derivada_central_2 = zeros(size(h_values_2));
    error_central_2 = zeros(size(h_values_2));

    fprintf('h\t\tDerivada Central\tValor Exacto\t\tError Absoluto\tlog(h)\t\tlog(Error)\n');
    fprintf('---\t\t----------------\t-------------\t\t---------------\t------\t\t----------\n');

    for i = 1 : length(h_values_2)
        h = h_values_2(i);

        % Calcular segunda derivada central
        deriv_central_2(i) = (f(x0-h) - 2*f(x0) + f(x0 + h)) / (h^2);

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

    % GENERACIÓN DE GRÁFICOS DE CONVERGENCIA 
    fprintf('\n=== GRÁFICOS DE CONVERGENCIA ===\n');

    figure('Position', [100, 100, 1200, 800]);

    % Gráfico 1: Diferencias Adelantadas
    subplot(2, 2, 1);
    loglog(h_values, error_derivada, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    loglog(h_values, h_values, 'r--', 'LineWidth', 1.5);  % Línea de referencia O(h)
    xlabel('h (Δx)');
    ylabel('Error Absoluto');
    title('Convergencia - Diferencias Adelantadas');
    legend('Error Adelantado', 'O(h)', 'Location', 'northeast');
    grid on;

    % Gráfico 2: Diferencias Centrales (Primera Derivada)
    subplot(2, 2, 2);
    loglog(h_values, error_central_1, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    loglog(h_values, h_values.^2, 'r--', 'LineWidth', 1.5);  % Línea de referencia O(h²)
    xlabel('h (Δx)');
    ylabel('Error Absoluto');
    title('Convergencia - Diferencias Centrales (1ª Derivada)');
    legend('Error Central', 'O(h²)', 'Location', 'northeast');
    grid on;

    % Gráfico 3: Diferencias Centrales (Segunda Derivada)
    subplot(2, 2, 3);
    loglog(h_values_2, error_central_2, 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    loglog(h_values_2, h_values_2.^2, 'r--', 'LineWidth', 1.5);  % Línea de referencia O(h²)
    xlabel('h (Δx)');
    ylabel('Error Absoluto');
    title('Convergencia - Diferencias Centrales (2ª Derivada)');
    legend('Error Central', 'O(h²)', 'Location', 'northeast');
    grid on;

    % Gráfico 4: Comparación de Métodos
    subplot(2, 2, 4);
    loglog(h_values, error_derivada, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    loglog(h_values, error_central_1, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
    loglog(h_values_2, error_central_2, 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('h (Δx)');
    ylabel('Error Absoluto');
    title('Comparación de Métodos');
    legend('Diferencias Adelantadas', 'Diferencias Centrales (1ª Derivada)', 'Diferencias Centrales (2ª Derivada)', 'Location', 'northeast');
    grid on;
endfunction

function ej2 
    fprintf('=== EJERCICIO 2: DERIVACIÓN NUMÉRICA ===\n');
    x0 = 0;
    xf = 2*pi;
    f = @(x) cos(x); 
    df_exact = @(x) -sin(x);
    d2f_exact = @(x) -cos(x);

    n_values = [5, 10, 20];

    % Vectores para análisis de convergencia
    delta_x = zeros(1, length(n_values));
    errores_primera = zeros(1, length(n_values));
    errores_segunda = zeros(1, length(n_values));
    fprintf('N\t\th (Δx)\t\tError 1ª Derivada\tError 2ª Derivada\n');
    fprintf('---\t\t-------\t\t----------------\t\t----------------\n');
    for i = 1 : length(n_values)
        n = n_values(i);
        h = (xf - x0) / n;
        delta_x(i) = h;

        % Discretización del intervalo
        x = linspace(x0, xf, n+1);
        y = f(x);

        % Calcular derivada  
        df_numerica = zeros(size(x));
        d2f_numerica = zeros(size(x));

        % Derivada primera
        for j = 2 : n
            df_numerica(j) = (y(j+1) - y(j-1)) / (2*h);
        endfor
        % Derivada segunda
        for j = 2 : n
            d2f_numerica(j) = (y(j+1) - 2*y(j) + y(j-1)) / (h^2);
        endfor

        % Calcular errores comparando con valores exactos
        df_numerica = df_exact(x);
        d2f_numerica = d2f_exact(x);

        %Error máximo en derivada primera (solo puntos intermedios)
        errores_primera(i) = max(abs(df_numerica(2:n) - df_exact(2:n)));
        %Error máximo en derivada segunda (solo puntos intermedios)
        errores_segunda(i) = max(abs(d2f_numerica(2:n) - d2f_exact(2:n)));
        
        % Mostrar resultados
        fprintf('%d\t\t%.5f\t\t%.2e\t\t%.2e\n', n, h, errores_primera(i), errores_segunda(i));

        % Crear gráficos para cada valor de N
        if i == 1  % Solo crear figura para el primer caso (N=5) como ejemplo
            figure('Position', [100, 100, 1200, 800]);
            
            % Gráfico 1: Función original
            subplot(2,2,1);
            plot(x, y, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            x_fine = linspace(x0, xf, 1000);
            plot(x_fine, f(x_fine), 'r-', 'LineWidth', 1.5);
            xlabel('x');
            ylabel('y = cos(x)');
            title(sprintf('Función Original (N=%d)', n));
            legend('Discreta', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 2: Primera derivada
            subplot(2,2,2);
            plot(x(2:n), df_numerica(2:n), 'go-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(x_fine, df_exact(x_fine), 'r-', 'LineWidth', 1.5);
            xlabel('x');
            ylabel("y' = -sin(x)");
            title(sprintf('Primera Derivada (N=%d)', n));
            legend('Numérica', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 3: Segunda derivada
            subplot(2,2,3);
            plot(x(2:n), d2f_numerica(2:n), 'mo-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(x_fine, d2f_exact(x_fine), 'r-', 'LineWidth', 1.5);
            xlabel('x');
            ylabel("y'' = -cos(x)");
            title(sprintf('Segunda Derivada (N=%d)', n));
            legend('Numérica', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 4: Comparación de errores
            subplot(2,2,4);
            loglog(delta_x(1:i), errores_primera(1:i), 'go-', 'LineWidth', 2, 'MarkerSize', 8);
            hold on;
            loglog(delta_x(1:i), errores_segunda(1:i), 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
            loglog(delta_x(1:i), delta_x(1:i).^2, 'r--', 'LineWidth', 1.5);
            xlabel('h (Δx)');
            ylabel('Error Máximo');
            title('Convergencia de Errores');
            legend('1ª Derivada', '2ª Derivada', 'O(h²)', 'Location', 'northeast');
            grid on;
            
            title(sprintf('Derivación Numérica de cos(x) - N=%d', n));
        endif
    endfor
endfunction

function ej3 
    fprintf('=== EJERCICIO 3: DERIVACIÓN NUMÉRICA ===\n');
    x0 = 0;
    xf = 2;
    f = @(t) (1/2) * (1 - exp(-2*t)); 
    df_exact = @(t) exp(-2*t);
    d2f_exact = @(t) -2*exp(-2*t);

    n_values = [5, 10, 20];

    % Vectores para análisis de convergencia
    delta_x = zeros(1, length(n_values));
    errores_primera = zeros(1, length(n_values));
    errores_segunda = zeros(1, length(n_values));
    fprintf('N\t\th (Δx)\t\tError 1ª Derivada\tError 2ª Derivada\n');
    fprintf('---\t\t-------\t\t----------------\t\t----------------\n');
    for i = 1 : length(n_values)
        n = n_values(i);
        h = (xf - x0) / n;
        delta_x(i) = h;

        % Discretización del intervalo
        t = linspace(x0, xf, n+1);
        y = f(t);

        % Calcular derivada  
        df_numerica = zeros(size(t));
        d2f_numerica = zeros(size(t));

        % Derivada primera
        for j = 2 : n
            df_numerica(j) = (y(j+1) - y(j-1)) / (2*h);
        endfor
        % Derivada segunda
        for j = 2 : n
            d2f_numerica(j) = (y(j+1) - 2*y(j) + y(j-1)) / (h^2);
        endfor

        % Calcular errores comparando con valores exactos
        df_numerica = df_exact(t);
        d2f_numerica = d2f_exact(t);

        %Error máximo en derivada primera (solo puntos intermedios)
        errores_primera(i) = max(abs(df_numerica(2:n) - df_exact(2:n)));
        %Error máximo en derivada segunda (solo puntos intermedios)
        errores_segunda(i) = max(abs(d2f_numerica(2:n) - d2f_exact(2:n)));
        
        % Mostrar resultados
        fprintf('%d\t\t%.5f\t\t%.2e\t\t%.2e\n', n, h, errores_primera(i), errores_segunda(i));

        % Crear gráficos para cada valor de N
        if i == 1  % Solo crear figura para el primer caso (N=5) como ejemplo
            figure('Position', [100, 100, 1200, 800]);
            
            % Gráfico 1: Función original
            subplot(2,2,1);
            plot(t, y, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            t_fine = linspace(x0, xf, 1000);
            plot(t_fine, f(t_fine), 'r-', 'LineWidth', 1.5);
            xlabel('t');
            ylabel('y = (1/2)(1 - e^{-2t})');  % CORREGIDO: etiqueta correcta
            title(sprintf('Función Original (N=%d)', n));
            legend('Discreta', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 2: Primera derivada
            subplot(2,2,2);
            plot(t(2:n), df_numerica(2:n), 'go-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(t_fine, df_exact(t_fine), 'r-', 'LineWidth', 1.5);
            xlabel('t');
            ylabel("y' = e^{-2t}");  % CORREGIDO: etiqueta correcta
            title(sprintf('Primera Derivada (N=%d)', n));
            legend('Numérica', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 3: Segunda derivada
            subplot(2,2,3);
            plot(t(2:n), d2f_numerica(2:n), 'mo-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(t_fine, d2f_exact(t_fine), 'r-', 'LineWidth', 1.5);
            xlabel('t');
            ylabel("y'' = -2e^{-2t}");  % CORREGIDO: etiqueta correcta
            title(sprintf('Segunda Derivada (N=%d)', n));
            legend('Numérica', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 4: Comparación de errores
            subplot(2,2,4);
            loglog(delta_x(1:i), errores_primera(1:i), 'go-', 'LineWidth', 2, 'MarkerSize', 8);
            hold on;
            loglog(delta_x(1:i), errores_segunda(1:i), 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
            loglog(delta_x(1:i), delta_x(1:i).^2, 'r--', 'LineWidth', 1.5);
            xlabel('h (Δt)');
            ylabel('Error Máximo');
            title('Convergencia de Errores');
            legend('1ª Derivada', '2ª Derivada', 'O(h²)', 'Location', 'northeast');
            grid on;
            
            title(sprintf('Derivación Numérica de (1/2)(1-e^{-2t}) - N=%d', n));
        endif
    endfor
endfunction

function ej4 
    fprintf('=== EJERCICIO 4: DERIVACIÓN NUMÉRICA ===\n');
    x0 = 0;
    xf = 8;
    f = @(t) exp(-0.5*t).*sin(4*t); 
    df_exact = @(t) exp(-0.5*t).*(4*cos(4*t) - 0.5*sin(4*t));
    d2f_exact = @(t) exp(-0.5*t).*(-16*sin(4*t) - 4*cos(4*t) + 0.25*sin(4*t));

    n_values = [10, 50, 100];

    % Vectores para análisis de convergencia
    delta_x = zeros(1, length(n_values));
    error_primera = zeros(1, length(n_values));
    error_segunda = zeros(1, length(n_values));
    fprintf('N\t\th (Δx)\t\tError 1ª Derivada\tError 2ª Derivada\n');
    fprintf('---\t\t-------\t\t----------------\t\t----------------\n');
    for i = 1 : length(n_values)
        n = n_values(i);
        h = (xf - x0) / n;
        delta_x(i) = h;

        % Discretización del intervalo
        t = linspace(x0, xf, n+1);
        y = f(t);

        % Calcular derivada  
        df_numerica = zeros(size(t));
        d2f_numerica = zeros(size(t));

        % Derivada primera
        for j = 2 : n
            df_numerica(j) = (y(j+1) - y(j-1)) / (2*h);
        endfor
        % Derivada segunda
        for j = 2 : n
            d2f_numerica(j) = (y(j+1) - 2*y(j) + y(j-1)) / (h^2);
        endfor

        % Calcular errores comparando con valores exactos
        df_exacta = df_exact(t);
        d2f_exacta = d2f_exact(t);

        %Error máximo en derivada primera (solo puntos intermedios)
        error_primera(i) = max(abs(df_numerica(2:n) - df_exacta(2:n)));

        %Error máximo en derivada segunda (solo puntos intermedios)
        error_segunda(i) = max(abs(d2f_numerica(2:n) - d2f_exacta(2:n)));
        
        % Mostrar resultados
        fprintf('%d\t\t%.5f\t\t%.2e\t\t%.2e\n', n, h, error_primera(i), error_segunda(i));

        % Crear gráficos para cada valor de N
        if i == 1  % Solo crear figura para el primer caso (N=5) como ejemplo
            figure('Position', [100, 100, 1200, 800]);
            
            % Gráfico 1: Función original
            subplot(2,2,1);
            plot(t, y, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            t_fine = linspace(x0, xf, 1000);
            plot(t_fine, f(t_fine), 'r-', 'LineWidth', 1.5);
            xlabel('t');
            ylabel('y = e^{-0.5t}sin(4t)');  % CORREGIDO: etiqueta correcta
            title(sprintf('Función Original (N=%d)', n));
            legend('Discreta', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 2: Primera derivada
            subplot(2,2,2);
            plot(t(2:n), df_numerica(2:n), 'go-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(t_fine, df_exact(t_fine), 'r-', 'LineWidth', 1.5);
            xlabel('t');
            ylabel("y' = e^{-0.5t}(4cos(4t) - 0.5sin(4t))");  % CORREGIDO: etiqueta correcta
            title(sprintf('Primera Derivada (N=%d)', n));
            legend('Numérica', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 3: Segunda derivada
            subplot(2,2,3);
            plot(t(2:n), d2f_numerica(2:n), 'mo-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
            plot(t_fine, d2f_exact(t_fine), 'r-', 'LineWidth', 1.5);
            xlabel('t');
            ylabel("y'' = e^{-0.5t}(-16sin(4t) - 4cos(4t) + 0.25sin(4t))");  % CORREGIDO: etiqueta correcta
            title(sprintf('Segunda Derivada (N=%d)', n));
            legend('Numérica', 'Exacta', 'Location', 'northeast');
            grid on;
            
            % Gráfico 4: Comparación de errores
            subplot(2,2,4);
            loglog(delta_x(1:i), error_primera(1:i), 'go-', 'LineWidth', 2, 'MarkerSize', 8);
            hold on;
            loglog(delta_x(1:i), error_segunda(1:i), 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
            loglog(delta_x(1:i), delta_x(1:i).^2, 'r--', 'LineWidth', 1.5);
            xlabel('h (Δt)');
            ylabel('Error Máximo');
            title('Convergencia de Errores');
            legend('1ª Derivada', '2ª Derivada', 'O(h²)', 'Location', 'northeast');
            grid on;
            
            title(sprintf('Derivación Numérica de e^{-0.5t}sin(4t) - N=%d', n));
        endif
    endfor
endfunction