function save_figure_safely(figHandle, fileName)
    try
        [parentFolder,~,~] = fileparts(fileName);
        if ~isempty(parentFolder) && ~exist(parentFolder,'dir')
            mkdir(parentFolder);
        end
        apply_unified_publication_style(figHandle);
        exportgraphics(figHandle, fileName, 'Resolution', 600);
    catch
        try
            saveas(figHandle, fileName);
        catch
            warning('Could not save figure: %s', fileName);
        end
    end
end
