function H = HSingle(z, m, g, l)
    H = z(2)^2/(2*m) - m*g/l * cos(z(1));
end