function c
    x = [0.5 1.5 2.5 3.5]'
    x2 = x.^2;
    y = [1.3 4.2 5.7 8.2]'
    v1 = [1 1 1 1]'

    fi2 = [v1 x x2]
    fi1 = [v1 x]

    a2 = (fi2' * fi2)^-1 * fi2' * y
    a1 = (fi1' * fi1)^-1 * fi1' * y

    r2 = y - fi2 * a2
    r1 = y - fi1 * a1


    # Aproximacion p
    p2 = fi2 * a2
    p1 = fi1 * a1

    normEr = norm(r1,2)
    normEr = norm(r2,2)
    errorProm1 = mean(r1)
    errorProm2 = mean(r2)
    yProm = mean(y)
    pProm = mean(p1)
    pProm2 = mean(p2)
    modr = dot(r1,r1)

    figure(1)
    plot(x,y,'-.',x,p2,'.-r',x,p1,'.-g')
    legend("Y","p2",'p1')
    grid on

endfunction