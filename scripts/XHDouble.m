function XH = XHDouble(z, m1, m2, l1, l2, g)

    % Define h1
    h1 = @(q, p) (p(1) * p(2) * sin(q(1) - q(2))) / ...
        (l1 * l2 * (m1 + m2 * sin(q(1) - q(2))^2));
    
    % Define h2
    h2 = @(q, p) sin(2*(q(1)-q(2)))*(m2 * l2^2 * p(1)^2 + (m1+m2) * l1^2 * p(2)^2 ...
        - 2 * m2 * l1 * l2 * p(1) * p(2) * cos(q(1) - q(2))) / ...
        (2 * l1^2 * l2^2 * (m1 + m2 * sin(q(1) - q(2))^2)^2);
    
    % Define H_q
    H_q = @(q, p) [
        (m1 + m2) * g * l1 * sin(q(1)) + h1(q, p) - h2(q, p);
        m2 * g * l2 * sin(q(2)) - h1(q, p) + h2(q, p)
    ];
    
    % Define H_p
    H_p = @(q, p) [
        (l2 * p(1) - l1 * p(2) * cos(q(1) - q(2))) / ...
        (l1^2 * l2 * (m1 + m2 * sin(q(1) - q(2))^2));
        (-m2 * l2 * p(1) * cos(q(1) - q(2)) + (m1 + m2) * l1 * p(2)) / ...
        (m2 * l1 * l2^2 * (m1 + m2 * sin(q(1) - q(2))^2))
    ];
    
    XH = [H_p(z(1:2),z(3:4));-H_q(z(1:2),z(3:4))];

end