function graf_ej6
    y0 = 1;
    t0 = 0;
    Dt = 0.1;
    NDt = 35;
    w=1;

    y = zeros(NDt,1);
    t = zeros(NDt,1);
    exacta = zeros(NDt,1);
    errorY = zeros(NDt,1);

    t(1) = t0;
    y(1) = y0;

    function edo = f_edo(y,t)
        edo = - t * y;
    end

    # Runge Kutta 2do Orden
    for j=1 : NDt - 1
        k1 = Dt * f_edo(y(j),t(j));

        yg = y(j) + k1 / (2*w);
        tg = t(j) + Dt / (2*w);
        
        k2 = Dt * f_edo(yg,tg);

        y(j+1) = y(j) + (1-w) * k1 + (w * k2);
        t(j+1) = t(j) + Dt;
    end

    # Exacta 
    for j = 1:NDt 
        exacta(j) = exp(-(t(j)^2)/2);
    end

    fprintf("Exacta y(0,2)= %f - Aprox = %f\n", exacta(3),y(3));
    figure(1)
    plot(t,exacta,'-b',t,y,'.r')
    grid on
    legend("Exacta","Runge Kutta")
endfunction