function a

    # VECTOR Y (Traspuesto para que quede columna)
    y = [2.3 4.5 7.3];
    y = y';

    # Vectores fi1 y f2, traspuestos para que sean columnas
    fi1 = [1 2 3];

    fi2 = [1 2.5 4];

    fi = [fi1' fi2'];

    # Vector a 
    a = (fi' * fi)^-1 * fi' * y

    # Aproximacion p
    p = fi * a

    # Residuo r
    r = y - p
    modr = dot(r,r)

endfunction