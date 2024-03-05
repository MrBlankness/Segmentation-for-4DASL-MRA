% 重写一个由线性索引查找下标的函数，适合多维，输出数组代表下标
% matlab 自带输出结果必须指定下标个数，此处无需指定
function subarray = myind2sub(datasize, ind)
% input：     datasize    数据尺寸
%             ind         数据线性索引 整数
% output:     subarray    数组形式的下标索引
% suozi   2016.05.17 HIT
% 379786867  buaasuozi@126.com
% ind 判断
if ind ~= fix(ind)
    disp('输入的索引必须为整数')
    return
end

subarray = zeros(size(datasize));
rest = ind;
i=length(datasize);
while i > 0
    if i ~= 1
       tmpdivide = rest/prod(datasize(1:i-1));
        if tmpdivide == fix(tmpdivide) % 余数为0 
            subarray(i) = tmpdivide;
        else
            subarray(i) = floor(tmpdivide) + 1;
        end
        tmprest = rest - (subarray(i) - 1)*prod(datasize(1:i-1));
        if tmprest ~= 0
            rest = tmprest;
            % else   rest = rest;
        end
    else
        subarray(i) = rest;
    end
    i = i - 1;
end

end
