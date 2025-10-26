function graf_ej7
    y10=5;
    y20=3;
    t0 = 0;
    Dt = 0.02;
    NDt = 100;
    w=1;

    t = zeros(1,NDt);
    y = zeros(2,NDt);
    yMod = zeros(2,NDt);
    ya = zeros(2,1);
    k1 = zeros(2,1);
    exacta = zeros(1,NDt);
    errorAd = zeros(1,NDt);
    errorMod = zeros(1,NDt);

    t(1) = t0;
    y(1,1) = y10;
    y(2,1) = y20;
    yMod(1,1) = y10;
    yMod(2,1) = y20;

    function [edo] = f_edo(y,t)
        edo(1,1) = -10 * y(1) + 4 * y(2);
        edo(2,1) = -4 * y(1);
    end

    # Euler adelante
    for j=1:NDt -1
        ya = y(:,j);

        k1 = Dt * f_edo(ya,t(j));
        y(:,j+1) = ya + k1;
        t(j+1) = t(j) + Dt;
    end

    # Euler modificado
    for j=1:NDt - 1
        ya = yMod(:,j);

        k1 = Dt * f_edo(ya,t(j));

        yg = ya + k1/(2*w);
        tg = t(j) + Dt/(2*w);

        k2 = Dt * f_edo(yg,tg);

        yMod(:,j+1) = ya + (1-w)*k1 + w*k2;
        t(j+1) = t(j) + Dt;
    end

    # Exacta
    for i=1:NDt
        exacta(i) = (1/3) * exp(-2 * t(i)) + (14/3) * exp(-8 * t(i));
    end

    # Errores
    for j=1:NDt
        errorAd(j) = abs(exacta(j)-y(1,j));
        errorMod(j) = abs(exacta(j) - yMod(1,j));
    end

    figure(1)
    plot(t,y(1,:),'--r',t,exacta,'-b',t,yMod(1,:),'-g')
    legend("Euler adelante","Exacta","Euler modificado")
    grid on

    figure(2)
    plot(t,errorAd,'-r',t,errorMod,'-b')
    legend("Error Euler Adelante","Error Euler Modificado")
    grid on

    
endfunction