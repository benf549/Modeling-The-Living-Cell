function dydt = predatorprey(t, y, k1, k2, k4)
    dydt = zeros(2,1);
    dydt(1) = k1*y(1) - k2*y(1)*y(2);
    dydt(2) = k2*y(1)*y(2) - k4*y(2);
end
