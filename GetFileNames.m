%gets fileNames filtered by 'fileExtension' from 'directory' recursively
%usage: fNames = GetFileNames('D:\bla\', '.m');
%files=GetFileNames(cd,'.txt');
function [fileNames, pathNames]= GetFileNames(directory, fileExtension, varargin)
    
    contents    = dir(directory);
    directories = find([contents.isdir]);

    fileIndicies = ~[contents.isdir];    
    fileStructures = contents(fileIndicies);

    %get files
    fileNames = {}; pathNames = {};
    for i = 1 : 1 : length(fileStructures)
        %fileName = fullfile(directory, fileStructures(i).name); %带文件位置的文件名
        fileName = fullfile(fileStructures(i).name);
        %**************filter****************
        %[folder, name, extension] = fileparts(fileName);
        [~, ~, extension] = fileparts(fileName);
        if( strcmp(extension, fileExtension) )
            fileNames = cat(1, fileNames, fileName);
            pathNames = cat(1, pathNames, directory);
        end
    end

    %recurse down (directory tree)搜索子文件夹
    if strcmp(varargin, 'recurse')
        for idxDir = directories
            
            subDirectory  = contents(idxDir).name;
            fullDirectory = fullfile(directory, subDirectory);
            
            % ignore '.' and '..'
            if (strcmp(subDirectory, '.') || strcmp(subDirectory, '..'))
                continue;
            end
            
            % Recurse down
            fileNames = cat(1, fileNames, GetFileNames(fullDirectory, fileExtension));
            [~,temppath]=GetFileNames(fullDirectory, fileExtension);
            pathNames = cat(1, pathNames, temppath);
        end
    end
end