% ��дһ�����������������±�ĺ������ʺ϶�ά�������������±�
% matlab �Դ�����������ָ���±�������˴�����ָ��
function subarray = myind2sub(datasize, ind)
% input��     datasize    ���ݳߴ�
%             ind         ������������ ����
% output:     subarray    ������ʽ���±�����
% suozi   2016.05.17 HIT
% 379786867  buaasuozi@126.com
% ind �ж�
if ind ~= fix(ind)
    disp('�������������Ϊ����')
    return
end

subarray = zeros(size(datasize));
rest = ind;
i=length(datasize);
while i > 0
    if i ~= 1
       tmpdivide = rest/prod(datasize(1:i-1));
        if tmpdivide == fix(tmpdivide) % ����Ϊ0 
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
