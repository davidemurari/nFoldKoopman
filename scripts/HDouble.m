function H = HDouble(z, m1, m2, l1, l2, g)
    q = z(1:2);
    p = z(3:4);
    H = (m2 * l2^2 * p(1)^2 + (m1 + m2) * l1^2 * p(2)^2 - 2 * m2 * l1 * l2 * p(1) * p(2) * cos(q(1) - q(2))) ...
    / (2 * m2 * l1^2 * l2^2 * (m1 + m2 * sin(q(1) - q(2))^2)) ...
    - (m1 + m2) * g * l1 * cos(q(1)) - m2 * g * l2 * cos(q(2));
end
