%input vector
x = 0.1:1/22:1;
x_size = length(x);
y1 = (1 + 0.6*sin(2*pi*x/0.7)) + (0.3*sin(2*pi*x))/2;

% plot graph
plot(x,y1);

c1 = 0.2;
c2 = 0.9;
r1 = 0.15;
r2 = 0.2;
w1 = randn(1);
w2 = randn(1);
b = randn(1);
nu = 0.15;


e = 0;
for i = 1:1000
    E = zeros(1, x_size);
    for j = 1:x_size
        
        g1 = (exp(-(x(j)-c1)^2 / (2*r1^2)));
        g2 = (exp(-(x(j)-c2)^2 / (2*r2^2)));
        f = g1*w1 + g2*w2 + b;  

        y = f;
        T = y1(j);

        E(j) = T - y;

        w1 = w1 + nu*E(j)*g1;
        w2 = w2 + nu*E(j)*g2;
        b = b + nu*e;
    end
    
end

hold on;

y_test = zeros(1, x_size);
for i = 1:x_size
    g1 = exp(-(x(i)-c1)^2/(2*r1^2));
    g2 = exp(-(x(i)-c2)^2/(2*r2^2));

    y_test(i) = g1*w1 + g2*w2 + b; 
end

disp(y_test);
plot(x, y_test);
