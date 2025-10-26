function regresionlineal 
    phi1 = [1; 1; 1; 1]
    phi2 = [0.5; 1.5; 2.5; 3.5] 
    y = [1.3; 4.2; 5.7; 8.2]

    PHI = [phi1  phi2]
    #´
    a = inv(PHI'* PHI) * (PHI') * y
    #´
    p = PHI * a
    r = y - p
    modulo = r'*r

    % Discretización de la función y para graficar
    t0 = 0;
    tf = 5;
    dt = 0.1; % Distancia entre cada intervalo
    N = floor((tf- t0)/dt); % Cantidad de intervalos 
    t(1) = t0;
    u(1) =  a(1, 1) + a(2, 1)*t(1); 

    for i=1 : N + 1 
      t(i+1) = t(i) + dt;
      u(i+1) = a(1, 1) + a(2, 1)*t(i+1); 
    endfor
    figure(1);
    hold on; 
    plot(t, u, '-r', 'DisplayName', 'Ajuste Lineal');
    plot(phi2, y, 'ob', 'DisplayName', 'Datos');
    grid on;
    legend('show')

    % Mantener las figuras abiertas
    fprintf('\nPresiona Enter para cerrar las figuras...\n');
    pause;

endfunction

% Ejecutar la función
regresionlineal;
