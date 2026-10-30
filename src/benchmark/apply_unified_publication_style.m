function apply_unified_publication_style(figHandle)
    % Unified post-processing style for all publication figures.
    % Required format:
    %   - Times New Roman
    %   - grid off, box on
    %   - line thickness = 2
    %   - bold x/y labels and bold tick labels, font size 14
    %   - alphabet-only subplot labels positioned above each axes

    if nargin < 1 || isempty(figHandle)
        figHandle = gcf;
    end

    set(figHandle,'Color','w');

    axList = findall(figHandle,'Type','axes');
    for ia = 1:numel(axList)
        ax = axList(ia);
        try
            set(ax, 'FontName','Times New Roman', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'LineWidth',2, ...
                'Box','on');
            grid(ax,'off');

            if ~isempty(ax.XLabel)
                set(ax.XLabel, 'FontName','Times New Roman', ...
                    'FontSize',14, 'FontWeight','bold');
            end
            if ~isempty(ax.YLabel)
                set(ax.YLabel, 'FontName','Times New Roman', ...
                    'FontSize',14, 'FontWeight','bold');
            end
            if ~isempty(ax.ZLabel)
                set(ax.ZLabel, 'FontName','Times New Roman', ...
                    'FontSize',14, 'FontWeight','bold');
            end
            if ~isempty(ax.Title) && strlength(string(ax.Title.String)) > 0
                set(ax.Title, 'FontName','Times New Roman', ...
                    'FontSize',18, 'FontWeight','bold');
                ax.Title.Units = 'normalized';
                ax.Title.Position(2) = 1.08;
            end
        catch
            % Continue styling the remaining axes even if one object type fails.
        end
    end

    % Make all plotted curves and markers visually clear.
    lineList = findall(figHandle,'Type','Line');
    for il = 1:numel(lineList)
        try
            set(lineList(il),'LineWidth',2);
        catch
        end
    end

    % Make bar/histogram edges clear where supported.
    graphicTypes = {'Bar','Histogram','ConstantLine'};
    for it = 1:numel(graphicTypes)
        objList = findall(figHandle,'Type',graphicTypes{it});
        for io = 1:numel(objList)
            try
                set(objList(io),'LineWidth',2);
            catch
            end
        end
    end

    % Style legends and colorbar labels.
    legList = findall(figHandle,'Type','Legend');
    for il = 1:numel(legList)
        try
            set(legList(il),'FontName','Times New Roman', ...
                'FontSize',14,'FontWeight','bold','Box','on');
        catch
        end
    end

    cbList = findall(figHandle,'Type','ColorBar');
    for ic = 1:numel(cbList)
        try
            set(cbList(ic),'FontName','Times New Roman', ...
                'FontSize',14,'FontWeight','bold','LineWidth',2);
            cbList(ic).Label.FontName = 'Times New Roman';
            cbList(ic).Label.FontSize = 14;
            cbList(ic).Label.FontWeight = 'bold';
        catch
        end
    end

    txtList = findall(figHandle,'Type','Text');
    for itxt = 1:numel(txtList)
        try
            set(txtList(itxt),'FontName','Times New Roman','FontWeight','bold');
            if txtList(itxt).FontSize < 14
                txtList(itxt).FontSize = 14;
            end
        catch
        end
    end
end
