% Optimal Attitude Determination Algorithm Using Fast Singular Value Decomposition
% author: liuzhuohua@bupt.edu.cn; jin_wu_uestc@hotmail.com


function [UU, VV] = svd3_fast(A)

    A11 = A(1, 1);            A12 = A(1, 2);                A13 = A(1, 3);
    A21 = A(2, 1);            A22 = A(2, 2);                A23 = A(2, 3);
    A31 = A(3, 1);            A32 = A(3, 2);                A33 = A(3, 3);
                              
    B11 = A11 * A11 + A21 * A21 + A31 * A31;  
    B12 = A11 * A12 + A21 * A22 + A31 * A32; 
    B13 = A11 * A13 + A21 * A23 + A31 * A33;
    B22 = A12 * A12 + A22 * A22 + A32 * A32; 
    B23 = A12 * A13 + A22 * A23 + A32 * A33;
    B33 = A13 * A13 + A23 * A23 + A33 * A33;
    
    B12B23 = B12 * B23;  
    B13B22 = B13 * B22; 
    B13B12 = B13 * B12; 
    B11B23 = B11 * B23; 
    B12B12 = B12 * B12;
    B11B22 = B11 * B22;
    
    a = -(B11 + B22 + B33);
    b =   B11B22 + B11 * B33 + B22 * B33 - B23^2 - B12B12 - B13^2;
    c = -(B11B22 * B33 - B11B23 * B23 - B12B12 * B33 + 2 * B12B23 * B13 - B13 * B13B22);
    a3 = a / 3;
    
    p = b - a * a3; 
    q = (2 * a^3 - 9 * a * b + 27 * c) / 27;   
    p2 = p * p;

    T = -q * sqrt(-27 * p) / 2 / p2;
    theta3 = acos(T) / 3;
    pi23 = 2.094395102393195;
    
    x0 = 2 * sqrt(-p / 3); 
    x1 = x0 * cos(theta3); 
    x2 = x0 * cos(theta3 - pi23); 
    x3 = x0 * cos(theta3 + pi23);

    rs1 = x1 - a3; 
    rs2 = x2 - a3; 
    rs3 = x3 - a3;
    
    if rs1 < rs2, tmp = rs2; rs2 = rs1; rs1 = tmp; end
    if rs1 < rs3, tmp = rs3; rs3 = rs1; rs1 = tmp; end
    if rs2 < rs3, tmp = rs2; rs2 = rs3; rs3 = tmp; end
    
    
    v1 = [ B12B23 - B13B22 + rs1 * B13
           B13B12 - B11B23 + rs1 * B23
           (rs1 - B11) * (rs1 - B22) - B12B12 ];
    v2 = [ B12B23 - B13B22 + rs2 * B13
           B13B12 - B11B23 + rs2 * B23
           (rs2 - B11) * (rs2 - B22) - B12B12 ];
           
    v1 = v1 / norm(v1); 
    v2 = v2 / norm(v2);
    
    v3 = [ v1(2) * v2(3) - v1(3) * v2(2);
           v1(3) * v2(1) - v1(1) * v2(3);
           v1(1) * v2(2) - v1(2) * v2(1) ];
       
    d1 = sqrt(rs1); 
    d2 = sqrt(rs2); 
    d3 = sqrt(rs3);
    
    u1 = A * v1 / d1; 
    u2 = A * v2 / d2; 
    u3 = A * v3 / d3;
    
    UU = [u1, u2, u3];
    VV = [v1, v2, v3];
end
