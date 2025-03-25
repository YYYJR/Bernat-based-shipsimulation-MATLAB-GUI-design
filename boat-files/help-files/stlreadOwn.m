function varargout = stlreadOwn(file)
%`stlreadOwn` 函数通过读取 STL 文件的内容，并调用 `stlbinary` 函数解析几何数据，支持返回顶点、面信息和法向量。
% 它适用于从 STL 文件中导入几何数据到 MATLAB，方便后续的可视化或分析。

% STLREAD imports geometry from an STL file into MATLAB.
%    FV = STLREAD(FILENAME) imports triangular faces from the ASCII or binary
%    STL file idicated by FILENAME, and returns the patch struct FV, with fields
%    'faces' and 'vertices'.
%
%    [F,V] = STLREAD(FILENAME) returns the faces F and vertices V separately.
%
%    [F,V,N] = STLREAD(FILENAME) also returns the face normal vectors.
%
%    The faces and vertices are arranged in the format used by the PATCH plot
%    object.
% Copyright 2011 The MathWorks, Inc.
    if ~exist(file,'file')
        error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
               ', be sure to specify the full path to the file.'], file);
    end
    
    fid = fopen(file,'r');    
    if ~isempty(ferror(fid))
        error(lasterror); %#ok
    end
    
    M = fread(fid,inf,'uint8=>uint8');
    fclose(fid);
    
    [f,v,n] = stlbinary(M);
    
    varargout = cell(1,nargout);
    switch nargout        
        case 2
            varargout{1} = f;
            varargout{2} = v;
        case 3
            varargout{1} = f;
            varargout{2} = v;
            varargout{3} = n;
        otherwise
            varargout{1} = struct('faces',f,'vertices',v);
    end
end


