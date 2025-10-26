function b

    y = [2.3 4.5 7.3];
    y = y';

    fi1 = [1 2 3];
    fi1 = fi1'

    a = (fi1' * fi1)^-1 * fi1' * y
    p = fi1 * a
    r = y - p
    modr = r'*r 
endfunction