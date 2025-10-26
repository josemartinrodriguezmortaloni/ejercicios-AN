function graf_ej8
    y10 = 10;
    y20 = 0;
    t0 = 0;
    Dt = 5.5E-3;
    w = 1/2;
    NDt = 250;

    y = zeros(2,NDt);
    t = zeros(1,NDt);
    k1 = zeros(2,1);
    k2 = zeros(2,1);

    y(1,1) = y10;
    y(2,1) = y20;
    t(1) = t0;

    function [edo] = f_edo(y,t)
        edo(1,1) = -0.5 * y(1) + 2 * y(2);
        edo(2,1) = -2 * y(1) - 0.5 * y(2);
    end

    for j=1:NDt -1
        ya = y(:,j);

        k1 = Dt * f_edo(ya,t(j));

        yg = ya + k1/(2*w);
        tg = t(j) + Dt/(2*w);

        k2 = Dt * f_edo(yg,tg);

        y(:,j+1) = ya + (1-w)*k1 + w*k2;
        t(j+1) = t(j) + Dt;
    end

    fo

    figure(1)
    plot(t,y(1,:),'.b')
    grid on

    figure(2)
    plot(y(1,:),y(2,:),'.b')
    grid on


endfunction
