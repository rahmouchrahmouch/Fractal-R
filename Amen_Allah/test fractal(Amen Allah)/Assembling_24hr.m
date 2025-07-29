load B11.mat; % Load 15min data
[rows, cols] = size(data);

% Ensure the number of rows is an integer
num_days = floor(rows / 96);
new = zeros(num_days, 1); % Initialize as a column vector

k = 0;
for i = 1:96:num_days*96
    k = k + 1;
    new(k) = sum(data(i:i+95, 1));
end

C11 = new / 0.0625;
save C11.mat

