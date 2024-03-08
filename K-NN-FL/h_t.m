function result = h_t(x_i_j, t, d)
    % Initialize auxiliary variable 'aux' to 1.
    aux = 1;

    % Loop over each dimension of the input vector 'x_i_j'.
    for i = 1:d
        % Multiply 'aux' by a scaled and shifted sine function of each element of 'x_i_j'.
        % The scaling factor is based on 't' and the constant 'pi/2'.
        % The sine function is scaled by sqrt(2).
        aux = sqrt(2) * sin(t * (pi / 2) * x_i_j(i)) * aux;
    end

    % Assign the final value of 'aux' to 'result'.
    result = aux;
end
