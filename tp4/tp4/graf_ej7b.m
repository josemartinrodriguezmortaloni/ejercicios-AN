function graf_ej7b
    y10=5;
    y20=3;
    t0=0;
    Dt = 0.01;
    NDt = 291;
    w = 0.5;

    y = zeros(2,NDt);
    yE = zeros(2,NDt);
    t = zeros(1,NDt);
    k1 = zeros(2,1);
    k2 = zeros(2,1);

    y(1,1) = y10;
    y(2,1) = y20;

    yE(1,1) = y10;
    yE(2,1) = y20;
    t(1) = t0;

    function [edo] = f_edo(y,t)
        edo(1,1) = y(2);
        edo(2,1) = -100 * y(1);
    end

    # Euler
    for j=1:NDt - 1
        ya = yE(:,j);
        k1 = Dt * f_edo(ya,t(j));
        yE(:,j+1) = ya + k1;
        t(j+1) = t(j) + Dt;
    end
    # Runge Kutta 2do Orden
    for j=1: NDt - 1
        ya = y(:,j);
        k1 = Dt * f_edo(ya,t(j));

        tg = t(j) + Dt/(2*w);
        yg = ya + k1/(2*w);

        k2 = Dt * f_edo(yg,tg);

        y(:,j+1) = ya + (1-w) * k1 + (w * k2);
        t(j+1) = t(j) + Dt;
    end

    figure(1)
    plot(t,yE(1,:),'.r',t,y(1,:),'.b')
    grid on

    figure(2)
    plot(y(1,:),y(2,:),'.b')
    grid on
endfunction