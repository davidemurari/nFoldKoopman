function XH = XHSingle(z, m, g, l)
    H_q = @(q, p)  m*g/l * sin(q);
    H_p = @(q, p) p/m;
    
    XH = [H_p(z(1),z(2));-H_q(z(1),z(2))];

end