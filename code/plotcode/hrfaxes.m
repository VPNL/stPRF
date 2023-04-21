
    function hrfaxes()
        
        xticklabels([xticks./10]);
        xlabel('Time (s)');
        ylabel('amplitude (a.u)')
        tch_set_axes;
        xlim([0 200]);
        ylim([-0.5 1.2])
        xticks([0    50   100   150   200])
        xticklabels([xticks./10]);
        yticks([-0.5 0 0.5 1 ])
        axis square;
    end