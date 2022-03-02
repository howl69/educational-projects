function GenerateTable(n,m)
    P = rand(n,m);
    save('table.mat', 'P')
    figure(1)
    image(P,'CDataMapping','scaled')
    colorbar
    disp(P)
end