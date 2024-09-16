function folders = GetFolders(directory)


filenames = {directory.name};
DirFlags = [directory.isdir] & ~strcmp(filenames, '.') & ~strcmp(filenames, '..');
folders = {directory(DirFlags).name};
