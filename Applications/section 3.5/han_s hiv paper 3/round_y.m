function new_y = round_y(y, epsilon)
    lower_counts = floor(y / epsilon);
    upper_counts = ceil(y / epsilon);
    lower_ints = lower_counts .* epsilon;
    upper_ints = upper_counts .* epsilon;
    logic = abs(lower_ints - y) < abs(upper_ints - y);
    new_y = logic .* lower_ints + (1 - logic) .* upper_ints;
end
 
