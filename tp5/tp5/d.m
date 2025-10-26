function d
    x = [4 6 8 10 12]'
    x2 = [1/4^2 1/6^2 1/8^2 1/10^2 1/12^2]'
    y = [9.37 9.64 9.74 9.78 9.81]'
    v1 = [1 1 1 1 1]'


    fi1 = [v1 x2]

    a1 = (fi1' * fi1)^-1 * fi1' * y

    r1 = y - fi1 * a1

    conv = a1 * [1 1/Inf^2]

    # Aproximacion p
    p1 = fi1 * a1

    figure(1)
    plot(x,p1,'.-r')