%input vector
x = 0.1:1/22:1;
x_size = length(x);
y1 = (1 + 0.6*sin(2*pi*x/0.7)) + (0.3*sin(2*pi*x))/2;

% plot graph
plot(x,y1);

c1 = 0.2;
c2 = 0.87;
c3 = 0.5;
r1 = 0.1;
r2 = 0.1;
r3 = 0.3;
w1 = randn(1);
w2 = randn(1);
w3 = randn(1);
b = randn(1);
nu = 0.03;
% Gauss function

e = 0;
for i = 1:100
    E = zeros(1, x_size);
    for j = 1:x_size
        
        g1 = (exp(-(x(j)-c1)^2 / (2*r1^2)));
        g2 = (exp(-(x(j)-c2)^2 / (2*r2^2)));
        g3 = exp(-(x(j)-c3)^2 /(2*r3^2));
        f = g1*w1 + g2*w2 + g3*w3 + b;  

        y = sign(f);
        T = y1(j);

        E(j) = T - y;
        e = e + E(j)^2;
%         
%         c1 = c1 + nu*(x(j) - c1);
%         c2 = c2 + nu*(x(j) - c2);
    end
    
    e = e/x_size;
    for j = 1:x_size
        g1 = exp(-(x(j)-c1)^2 / (2*r1^2));
        g2 = exp(-(x(j)-c2)^2 /(2*r2^2));
        
        w1 = w1 + nu*e*x(j);
        w2 = w2 + nu*e*x(j);
        w3 = w3 + nu*e*x(j);
        b = b + nu*e;
    end
end

% disp(w1);
% disp(w2);
% disp(b);

hold on;

y_test = zeros(1, x_size);
for i = 1:x_size
    g1 = exp(-(x(i)-c1)^2/(2*r1^2));
    g2 = exp(-(x(i)-c2)^2/(2*r2^2));

    y_test(i) = g1*w1 + g2*w2 + g3*w3 + b; 
end
% disp(y1);
disp(y_test);
plot(x, y_test);

function s = sign(x)
    if x > 0
        s = 1;
    else
        s = -1;
    end
end