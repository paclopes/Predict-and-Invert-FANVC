function plot_xy_p2(x,y)
    M = size(y,2);
    if M > 1
        sorted = sort(y,2);
        plot(x,sorted(:,ceil(0.50*M)), 'color', 'black');
        hold on;
        plot(x,sorted(:,ceil(0.01*M)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,ceil(0.25*M)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,ceil(0.75*M)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,ceil(0.99*M)), 'color', 'black', 'LineStyle',':');
        hold off;
    else
        plot(x,y, 'color', 'black')
    end

end