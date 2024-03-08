function result = xi_i(x, i,d,b_i_ts)
    aux = 0;
    for t = 1:5
        aux = aux + (b_i_ts(t,i)/(t^2))*h_t(x, t,d);
    end
    result = aux;
end
