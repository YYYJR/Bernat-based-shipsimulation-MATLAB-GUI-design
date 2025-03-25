function tf = isbinary(A)
% ISBINARY uses the first line of an STL file to identify its format.
    if isempty(A) || length(A) < 5
        error('MATLAB:stlread:incorrectFormat', ...
              'File does not appear to be an ASCII or binary STL file.');
    end    
    if strcmpi('solid',char(A(1:5)'))
        tf = false; % ASCII
    else
        tf = true;  % Binary
    end
end


%`isbinary` 函数通过检查 STL 文件的前 5 个字节来判断文件是 ASCII 格式还是二进制格式。
% 它基于 ASCII STL 文件以 `solid` 开头的特性，通过简单的字符串比较实现格式判断。
% 这种方法简单高效，适用于快速判断 STL 文件的格式。