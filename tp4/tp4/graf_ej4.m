function graf_ej4
    y0=1/10;
    t0=0;
    Dt = 0.01;
    w = 0.5;
    NDt = 50;

    y=zeros(NDt,1);
    t=zeros(NDt,1);
    exacta=zeros(NDt,1);
    errorY=zeros(NDt,1);

    y(1) = y0;
    t(1) = t0;

    function edo = f_edo(y,t)
        edo = exp(-2 * t) - 2 * y;
    end

    # Runge Kutta 2do Orden
    for j=1: NDt - 1
        k1 = Dt * f_edo(y(j),t(j));

        tg = t(j) + Dt/(2*w);
        yg = y(j) + k1/(2*w);

        k2 = Dt * f_edo(yg,tg);

        y(j+1) = y(j) + (1-w) * k1 + (w * k2);
        t(j+1) = t(j) + Dt;
    end

    # Exacta
    for j=1: NDt 
        exacta(j) = (1/10) * exp(-2 * t(j)) + t(j) * exp(-2 * t(j));
    end

    # Error
    for j=1:NDt -1
        errorY(j) = abs(exacta(j) - y(j));
    end

    # Norma infinita
    normaInf = norm(errorY,Inf)
    maxY = max(exacta) * 0.01;
    normaInf < maxY

    figure(1)
    plot(t,exacta,'-b',t,y,'-r')
    legend("Exacta","Runge Kutta")
    grid on

    figure(2)
    plot(t,errorY,'-')
    legend("Error")
    grid on

endfunction