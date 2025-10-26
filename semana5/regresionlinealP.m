function regresionlinealP 
    phi1 = [1; 1; 1; 1];
    phi2 = [0.5; 1.5; 2.5; 3.5]; 
    phi3 = phi2.^2;
    y = [1.3; 4.2; 5.7; 8.2];

    PHI = [phi1  phi2 phi3];
    #´
    a = inv(PHI'* PHI) * (PHI') * y;
    #´
    p = PHI * a;
    r = y - p;
    modulo = r'*r;

    % Discretización de las funciones para graficar
    t0 = 0;
    tf = 10;
    dt = 0.01; % Distancia entre cada intervalo
    N = floor((tf- t0)/dt); % Cantidad de intervalos
    t(1) = t0;
    u_lineal(1) = a(1, 1) + a(2, 1)*t(1);
    u_parabola(1) = a(1, 1) + a(2, 1)*t(1) + a(3, 1)*t(1)^2;

    for i=1 : N + 1
      t(i+1) = t(i) + dt;
      u_lineal(i+1) = a(1, 1) + a(2, 1)*t(i+1);
      u_parabola(i+1) = a(1, 1) + a(2, 1)*t(i+1) + a(3, 1)*t(i+1)^2;
    endfor
    for j=1 : length(a)
      fprintf('%.5f\n', a(j, 1))
    endfor

    figure(1);
    hold on;
    plot(t, u_lineal, '-r', 'DisplayName', 'Ajuste Lineal');
    plot(t, u_parabola, '-g', 'DisplayName', 'Ajuste Parabólico');
    plot(phi2, y, 'ob', 'DisplayName', 'Datos');
    grid on;
    xticks(0:0.5:10);  % Mostrar valores en X cada 0.5 para mejor legibilidad
    xlim([0 4]);       % Limitar el rango X a donde están los datos
    xlabel('X');
    ylabel('Y');
    title('Regresión Lineal vs Parabólica');
    legend('show')

    % Mantener las figuras abiertas
    fprintf('\nPresiona Enter para cerrar las figuras...\n');
    pause;

endfunction

% Ejecutar la función
regresionlinealP;
