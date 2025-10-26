function convergencia 
    phi1 = [1; 1; 1];
    N = [4; 6; 8];
    phi2 = (1 ./ N.^2);
    y = [9.37; 9.64; 9.74];

    PHI = [phi1  phi2]
    #´
    a = inv(PHI'* PHI) * (PHI') * y
    #´
    p = PHI * a
    r = y - p
    modulo = r'*r

    % Discretización de la función para graficar la convergencia
    t_min = 0;
    t_max = max(phi2) * 1.2;  % Un poco más que el máximo 1/N²
    dt = t_max / 100;         % 100 puntos para suavidad
    num_points = floor(t_max / dt);
    t_conv(1) = t_min;
    u_conv(1) = a(1, 1) + a(2, 1) * t_conv(1);

    for i = 1 : num_points
      t_conv(i+1) = t_conv(i) + dt;
      u_conv(i+1) = a(1, 1) + a(2, 1) * t_conv(i+1);
    endfor

    % Valor límite cuando N -> infinito (1/N² -> 0)
    valor_limite = a(1, 1);

    figure(1);
    hold on;
    plot(t_conv, u_conv, '-r', 'DisplayName', 'Ajuste Lineal');
    plot(phi2, y, 'ob', 'DisplayName', 'Datos');
    % Línea horizontal para el valor límite (compatible con Octave)
    x_range = [min(t_conv), max(t_conv)];
    plot(x_range, [valor_limite, valor_limite], '--k', 'DisplayName', sprintf('Límite = %.4f', valor_limite));
    grid on;
    xlabel('1/N²');
    ylabel('Valor convergente');
    title('Análisis de Convergencia');
    legend('show');

    % Mostrar información de convergencia
    fprintf('\n=== ANÁLISIS DE CONVERGENCIA ===\n');
    fprintf('Valor límite (N -> infinito): %.6f\n', valor_limite);
    fprintf('Pendiente: %.6f\n', a(2, 1));
    fprintf('Error residual: %.6f\n', modulo);

    % Mantener las figuras abiertas
    fprintf('\nPresiona Enter para cerrar las figuras...\n');
    pause;

endfunction

% Ejecutar la función
convergencia;
