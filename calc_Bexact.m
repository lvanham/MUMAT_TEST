clear;
%%
R = 5;
B0 = [0 0 1];
mu_r = 500;

chi_m = mu_r-1;

Bexact = @(r) (chi_m)/(1+chi_m/3)*(R/norm(r))^3*( dot(B0,r)*r/norm(r)^2 - B0/3 ) + B0;

%%
min1 = [500 0 0];
max1 = [2000 0.5*pi 2*pi];
num_points = [151 7 1];
npoints = prod(num_points);

n_temp = 1;
x = zeros(1,npoints);
y = zeros(1,npoints);
z = zeros(1,npoints);

for i = 1:num_points(1)
    r = min1(1) + (i-1)*(max1(1)-min1(1))/(num_points(1)-1);
    for j = 1:num_points(2)
        theta = min1(2) + (j-1)*(max1(2)-min1(2))/(num_points(2)-1);
        x(n_temp) = r*sin(theta);
        y(n_temp) = 0;
        z(n_temp) = r*cos(theta);
        n_temp = n_temp + 1;
    end
end

%%

B = zeros(length(x),3);
r = zeros(length(x),3);
for i = 1:length(x)
    r(i,:) = [x(i) y(i) z(i)]/100;
    B(i,:) = Bexact(r(i,:));
end

save Bexact.dat B -ascii