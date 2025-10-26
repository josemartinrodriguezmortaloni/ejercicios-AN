function graf_ej3
    y0 = 4;
    t0 = 0;
    Dt = 0.25;
    NDt = 40;
    w = 1/2;

    yEuler=zeros(NDt,1);
    yMej=zeros(NDt,1);
    yMod=zeros(NDt,1);
    t=zeros(NDt,1);
    exacta=zeros(NDt,1);

    yEuler(1) = y0;
    yMej(1) = y0;
    yMod(1) = y0;
    t(1) = t0;

    function edo = f_edo(y,t)
        edo = (t/2) - (1/2)*y;
    end

    # Euler
    for j=1:NDt -1
        k1 = Dt * f_edo(yEuler(j),t(j));

        yEuler(j+1) = yEuler(j) + k1;
        t(j+1) = t(j) + Dt;
    end

    # Euler Modificado
    for j=1:NDt - 1
        k1 = Dt * f_edo(yMod(j),t(j));

        yg = yMod(j) + k1/(2*w);
        tg = t(j) + Dt/(2*w);

        k2 = Dt * f_edo(yg,tg);

        yMod(j+1) = yMod(j) + ((1-w) * k1) + (w * k2);
        t(j+1) = t(j) + Dt;
    end

    # Euler Mejorado
    for j=1: NDt - 1
        k1 = Dt * f_edo(yMej(j),t(j));
        
        yg = yMej(j) + k1;
        tg = t(j) + Dt;

        k2 = Dt * f_edo(yg,tg);

        yMej(j+1) = yMej(j) + (k1+k2)/2;
        t(j+1) = t(j) + Dt;
    end

    # Exacta
    for j=1:NDt
        exacta(j) = 6*exp(-t(j)/2) - 2 + t(j);
    end

    # Errores
    errorEuler = abs(exacta(3)-yEuler(3))
    errorMod = abs(exacta(3)-yMod(3))
    errorMej = abs(exacta(3)-yMej(3))

    for j=1:NDt
        errorEuler(j) = abs(exacta(j)-yEuler(j));
        errorMod(j) = abs(exacta(j)-yMod(j));
        errorMej(j) = abs(exacta(j)-yMej(j));
    end

    figure(1)
    plot(t,exacta,'--',t,yEuler,'-r',t,yMod,"-b",t,yMej,'-.g')
    legend("Exacta","Euler","Modificado","Mejorado")
    grid on

    figure(2)
    plot(t,errorEuler,'-r',t,errorMod,"-b",t,errorMej,'-.g')
    legend("Euler","Modificado","Mejorado")
    grid on
endfunction