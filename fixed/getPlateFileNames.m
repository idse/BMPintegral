function [fileNames,wellNames] = getPlateFileNames(listing, conditions)

    fileNames = cell(size(conditions));
    wellNames = cell(size(conditions));

    for i = 1:numel(listing)

        fileName = listing(i).name;
        % ^ means starts with
        [~, matches] = strsplit(fileName,'(_| )[ABCDEFGH][0-9]+|^[ABCDEFGH][0-9]+',...
                                        'DelimiterType','RegularExpression');
        if ~isempty(matches)
            wellName = matches{1};
            if contains(wellName,'_') || contains(wellName,' ')
                wellName = wellName(2:end); % strip off the underscore or space
            end
            rowIdx = wellName(1) - 'A' + 1;
            colIdx = str2num(wellName(2:end));

            if isempty(fileNames{rowIdx,colIdx})
                fileNames{rowIdx,colIdx} = {fileName};
            else
                fileNames{rowIdx,colIdx} = [fileNames{rowIdx,colIdx}, fileName];
            end
            wellNames{rowIdx,colIdx} = wellName;
        end
    end

end