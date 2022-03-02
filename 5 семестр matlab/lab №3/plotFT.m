function res = plotFT(hFigure, fHandle, fFTHandle, step, inpLimVec, outLimVec)
    figure(hFigure);
    if isfield(hFigure.UserData, 'inpLimVec')
        a = hFigure.UserData.inpLimVec(1);
        b = hFigure.UserData.inpLimVec(2);
    else
        a = inpLimVec(1);
        b = inpLimVec(2);
    end
    if isfield(hFigure.UserData, 'outLimVec')
        c = hFigure.UserData.outLimVec(1);
        d = hFigure.UserData.outLimVec(2);
    elseif length(outLimVec) > 1 
        c = outLimVec(1);
        d = outLimVec(2);
    else
        c = -1/(2*step);
        d = 1/(2*step);
    end
    info = get(hFigure, 'UserData');
    if ~isempty(info) && isfield(hFigure.UserData, 'ax_handles')
        delete(info.ax_handles);
    else
        info = struct('first_axises', [], 'second_axises', []);
    end
    T = b - a;
    N = round(T / step);
    step = T / N;
    x = a:step:b;
    y = fHandle(x);
    y(isnan(y)) = 0;
    n = 0;
    if a > 0
        n = round(a/T);
    elseif b < 0
        n = -round(b/T);
    end
    border = n*T;
    ind = find(x <= border, 1, 'last' );
    Y(1:N+1-ind) = y(ind+1:N+1);
    Y(N+1-ind+1:N+1) = y(1:ind);
    Yft = fft(Y)*T/(N+1);
    stepft = 2*pi/T;
    Tft = stepft*(N+1);
    lbord = floor(c / Tft)*Tft;
    rbord = ceil(d / Tft)*Tft;
    counter = ceil(d / Tft) - floor(c / Tft);
    
    Xft = lbord:stepft:rbord - stepft;
    Yft = repmat(Yft(1:end), 1, counter);
    subplot(2, 1, 1);
    hold on;
    plot(Xft(1:end), real(Yft), 'b');
    legend('Апроксимация преобразования Фурье');
    if ~isempty(fFTHandle)
        res_f = fFTHandle(linspace(c, d, 1000));
        plot(linspace(c, d, 1000), real(res_f), 'r', 'DisplayName', 'Аналитически вычисленное преобразование Фурье');
        if isfield(hFigure.UserData, 'realXLim') 
            ylim([min(real(Yft))*1.2 max(real(Yft))*1.2])
        else
            axis([c d min([real(Yft) real(res_f)])   max([real(Yft) real(res_f)])]);
        end
    elseif isfield(hFigure.UserData, 'realXLim') 
        ylim([min(real(Yft))*1.2 max(real(Yft))*1.2])
    else
        axis([c d min(real(Yft)) max(real(Yft))]);
    end
    ylabel('Re F(\lambda)');
    xlabel('\lambda');
    title('Вещественная ось');
    info.first_axises = gca;  
    subplot(2, 1, 2);
    hold on;
    plot(Xft(1:end), imag(Yft), 'b');
    legend('Апроксимация преобразования Фурье');
    if ~isempty(fFTHandle)
        res_f = fFTHandle(linspace(c, d, 1000));
        plot(linspace(c, d, 1000), imag(res_f), 'r', 'DisplayName', 'Аналитически вычисленное преобразование Фурье');
        if isfield(hFigure.UserData, 'imagXLim') 
            ylim([min(imag(Yft))*1.2 max(imag(Yft))*1.2])
        else
            axis([c d min([imag(Yft) imag(res_f)])   max([imag(Yft) imag(res_f)])]);
        end
    elseif isfield(hFigure.UserData, 'imagXLim') 
        ylim([min(imag(Yft))*1.2 max(imag(Yft))*1.2])
    else
        axis([c d min(imag(Yft)) max(imag(Yft))]);
    end
    xlim([c d])
    title('Мнимая ось');
    ylabel('Im F(\lambda)');
    xlabel('\lambda');
    info.second_axises = gca;
    set(hFigure,'UserData',info);
    res.inpLimVec = [a b];
    res.outLimVec = [c d];
end