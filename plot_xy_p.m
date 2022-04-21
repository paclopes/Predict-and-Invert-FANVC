function plot_xy_p(x,y)
    M = size(y,2);
    if M > 1
        sorted = sort(y,2);

        plot(x,sorted(:,round(0.50*(M-1)+1)), 'color', 'black');
        hold on;
        plot(x,sorted(:,round(0.1*(M-1)+1)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,round(0.25*(M-1)+1)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,round(0.75*(M-1)+1)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,round(0.9*(M-1)+1)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,round(0.99*(M-1)+1)), 'color', [0.5 0.5 0.5], 'LineStyle',':');
        hold off;
    else
        plot(x,y, 'color', 'black')
    end

end