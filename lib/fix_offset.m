function fixed_data = fix_offset(data, signal)

edges = diff(signal);
rising = find(edges == 1);
falling = find(edges == -1);
num_rising = size(rising, 1);
num_falling = size(falling, 1);

assert(num_rising == num_falling, "edges not matching, failed to fix offset");

fixed_data = data;

for i = 1:num_rising
    fixed_data(rising(i):falling(i)) = fixed_data(rising(i):falling(i)) - data(rising(i));
end
end
