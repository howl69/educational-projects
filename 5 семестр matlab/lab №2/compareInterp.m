function compareInterp(x, xx, f)
    near_int = interp1(x,f(x),xx, 'nearest');
    lin_int = interp1(x,f(x),xx, 'linear');
    spl_int = interp1(x,f(x),xx, 'spline');
    pchip_int = interp1(x,f(x),xx, 'pchip');
    plot(xx, f(xx), xx, near_int, ':.', xx, lin_int, ':.', xx, spl_int, ':.', xx, pchip_int, ':.');
    legend('f(xx)', 'nearest', 'linear', 'spline', 'pchip')
    title('Interpolation methods');
end