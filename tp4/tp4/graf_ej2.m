function graf_ej2()
    t0 = 0;
    u0 = 2;
    w=1/2;

    Dt=0.25;  % incremento de tiempo
    NDt=8;  % cantidad de Dt a realizar

    % Dimensionamiento
    t=zeros(NDt,1);   % vector para tiempo
    t2 = zeros(NDt,1);
    u=zeros(NDt,1);   % vector para Euler Modificado
    y = zeros(NDt,1);  % vector para Euler Mejorado (corregido a columna)
    errorAbsMod=zeros(NDt,1);  % error absoluto Euler Modificado
    errorAbsMej=zeros(NDt,1);  % error absoluto Euler Mejorado
    exacta=zeros(NDt,1); % vector para la solución exacta

    u(1) = u0;
    t(1) = t0;
    t2(1) = t0;
    y(1) = u0;

    % Euler modificado
    for j=1:NDt - 1
        k1 = Dt * (2 * u(j) - 2 * t(j) - 1);

        ug = u(j) + k1/(2*w);
        tg = t(j) + Dt/(2*w);
        
        k2 = Dt * (2 * ug - 2 * tg - 1);

        t(j+1)  = t(j)+Dt;
        u(j+1)  = u(j) + (1-w)*k1 + w*k2;
    end

    % Euler Mejorado
    for j=1:NDt - 1
        k1 = Dt * (2 * y(j) - 2 * t2(j) - 1);
        yg = y(j) + k1;
        tg = t2(j) + Dt;
        k2 = Dt * (2 * yg - 2 * tg - 1);
        y(j+1) = y(j) + (k1 + k2)/2;
        t2(j+1) = t2(j) + Dt;
    end

    for j=1:NDt 
        exacta(j) = exp(2*t(j))+t(j)+1;
    end

    % Cálculo de errores absolutos
    for j=1:NDt
        errorAbsMod(j) = abs(u(j)-exacta(j));     % Error Euler Modificado
        errorAbsMej(j) = abs(y(j)-exacta(j));     % Error Euler Mejorado
    end

    fprintf("Euler modificado en f(0.5)= %f\n",u(3))
    fprintf("Euler mejorado en f(0.5)= %f\n",y(3))
    fprintf("Exacta en f(0.5)= %f\n",exacta(3))
    fprintf("Error absoluto Euler Modificado en t=0.5: %f\n", errorAbsMod(3))
    fprintf("Error absoluto Euler Mejorado en t=0.5: %f\n", errorAbsMej(3))

    % Gráfica de las soluciones
    figure(1)
    plot(t,u,"b.",t,exacta,"--g",t,y,'-r')
    title('Comparación de métodos numéricos vs solución exacta')
    xlabel('Tiempo (t)')
    ylabel('u(t)')
    legend("Euler modificado", "Exacta", "Euler mejorado")
    grid on
    
    % Gráfica de errores absolutos
    figure(2)
    plot(t,errorAbsMod,"r-", t,errorAbsMej,"b.");
    title('Error absoluto de los métodos numéricos')
    xlabel('Tiempo (t)')
    ylabel('Error absoluto')
    legend("Error Euler Modificado", "Error Euler Mejorado")
    grid on

endfunction
