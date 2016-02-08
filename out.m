Array = csvread('out.csv');
col1  = Array(:, 1);
col2  = Array(:, 2);
col3  = Array(:, 3);

scatter3(col1, col2, col3);
xlabel('nH');
ylabel('nH_p');
zlabel('dnH_p/dt');