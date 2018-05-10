function y = rv_ft(x, a, k, phi, m, b)
    rv = importdata('RV_tot.txt');
    y = a * sin(x/100*k+phi) * (1-m) + b + m * rv;
end 