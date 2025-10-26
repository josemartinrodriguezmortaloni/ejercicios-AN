function graf_ej5
    y0 = 1;
    t0 = 0;
    Dt = 0.01;
    NDt = 100;
    w=1;

    y = zeros(NDt,1);
    t = zeros(NDt,1);
    errorY = zeros(NDt,1);
    exacta = zeros(NDt,1);

    y(1) = y0;
    t(1) = t0;

    function edo = f_edo(y,t)
        edo = 2 * t * y^2;
    end

    for j=1:NDt -1
        k1 = Dt * f_edo(y(j),t(j));
        yg = y(j) + 1/(2*w) * k1;
        tg = t(j) + Dt/(2*w);

        k2 = Dt * f_edo(yg,tg);

        y(j+1) = y(j) + (1-w)* k1 + (w * k2);
        t(j+1) = t(j) + Dt;
    end

    # Exacta
    for k=1:NDt
        exacta(k) = 1/(1-t(k)^2);
    end

    # Error
    for k=1: NDt
        errorY(k) = abs(y(k) - exacta(k));
    end

    # Norma
    normaInf = norm(errorY,Inf)

    figure(1)
    plot(t,exacta,'-r',t,y,'--b')
    legend("exacta","runge kutta")
    grid on

    figure(2)
    plot(t,errorY,'-r')
    legend("ERROR")
    grid on


endfunction