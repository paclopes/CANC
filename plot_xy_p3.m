function plot_xy_p3(x,y)
    co = colororder;
    color1 = co(1,:);
    M = size(y,2);
    if M > 1
        sorted = sort(y,2);
        plot(x,sorted(:,ceil(0.50*M)), 'color', color1);
        hold on;
        plot(x,sorted(:,ceil(0.01*M)), 'LineStyle',':', 'color', color1);
        plot(x,sorted(:,ceil(0.25*M)), 'LineStyle',':', 'color', color1);
        plot(x,sorted(:,ceil(0.75*M)), 'LineStyle',':', 'color', color1);
        plot(x,sorted(:,ceil(0.99*M)), 'LineStyle',':', 'color', color1);
        hold off;
    else
        plot(x,y)
    end
end