function graf_ej9
    y10 = 0;
    y20 = 0;
    t0 = 0;
    Dt = 1E-2;
    NDt = 1500;
    w=1/2;

    y = zeros(2,NDt);
    t = zeros(1,NDt);
    k1 = zeros(2,1);
    k2 = zeros(2,1);

    function [edo] = f_edo(y,t)
        edo(1,1) = y(2) + sin(3*t) * 0;
        edo(2,1) = -4 * y(1) + sin(3*t) * 10;
    end

    for j=1:NDt - 1
        ya = y(:,j);

        k1 = Dt * f_edo(ya,t(j));

        yg = ya + k1/(2*w);
        tg = t(j) + Dt/(2*w);

        k2 = Dt * f_edo(yg,tg);

        y(:,j+1) = ya + (1-w) * k1 + w*k2;
        t(j+1) = t(j) + Dt;
    end
    
    figure(1)
    plot(t,y(1,:),'.b')
    grid on

    figure(2)
    plot(y(1,:),y(2,:),'.r')
    grid on

endfunction
