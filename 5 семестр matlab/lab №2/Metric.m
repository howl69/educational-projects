function metric = Metric(typeOfMetric, x, f, fn, n)
    if typeOfMetric
        metric = max(abs(f(x) - fn(x, n))); 
    else
        metric = (trapz(x, ((f(x) - fn(x, n)).^2)))^0.5;
    end
end