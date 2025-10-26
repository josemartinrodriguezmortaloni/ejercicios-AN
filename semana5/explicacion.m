function explicacion
    phi1 = [1; 1; 1; 1]
    phi2 = [0.5; 1.5; 2.5; 3.5] 
    y = [1.3; 4.2; 5.7; 8.2]

    PHI = [phi1  phi2]
    #´
    a = inv(PHI'* PHI) * (PHI') * y
    #´
    p = PHI * a
    r = y - p
    modulo = r'*r
endfunction

% Ejecutar la función
explicacion;
