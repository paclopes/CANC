function plot_xy_p4(x,y,i)
    co = colororder;
    color1 = co(1,:);
    M = size(y,2);
    if M > 1
        p = [0.01, 0.25, 0.5, 0.75, 0.99];
        sorted = sort(y,2);
        for j = 1:length(p)
            y = sorted(:, ceil(p(j)*M));
            if p(j) == 0.5
                ls = '-';
            else
                ls = ':';
            end
            plot(x, y, 'LineStyle', ls, 'color', color1);
            text(x(i), y(i), [num2str(100*p(j)), '\%'], 'VerticalAlignment', 'bottom');
            if j == 1
                hold on;
            end
            if j == length(p)
                hold off;
            end
        end
    else
        plot(x,y)
    end
end