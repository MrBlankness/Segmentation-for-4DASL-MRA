function sum_length = cal_skel_length(skel)
    addpath('./skel2graph3d-matlab-master');
    addpath('./skeleton3d-matlab-master');
    D = bwconncomp(skel, 26);
    disp('连通区域数目');
    disp(length(D.PixelIdxList));
    % volshow(skel);
    w = size(skel,1);
    l = size(skel,2);
    h = size(skel,3);
    [~,node,link] = Skel2Graph3D(skel,0);

    sum_length = 0;
    temp_x = 0;
    temp_y = 0;
    temp_z = 0;

%     figure
%     axis square
%     hold on
    for i = 1:length(link)
        x = zeros(length(link(i).point), 1);
        y = zeros(length(link(i).point), 1);
        z = zeros(length(link(i).point), 1);
        for j = 1 : length(link(i).point)
            if j ~= 1 && j ~= length(link(i).point)
                z(j) = ceil(link(i).point(j) / (w * l));
                y(j) = ceil((link(i).point(j) - ((z(j) - 1) * w * l)) / w);
                x(j) = link(i).point(j) - (z(j) - 1) * w * l - (y(j) - 1) * w;
            else
                judge = 0;
                for k = 1:length(node)
                    for kk = 1:length(node(k).idx)
                        if node(k).idx(kk) == link(i).point(j)
                            x(j) = node(k).comx;
                            y(j) = node(k).comy;
                            z(j) = node(k).comz;
                            judge = 1;
                            break
                        end
                    end
                    if judge == 1
                        break
                    end
                end
            end
        end
        xx = [];
        yy = [];
        zz = [];
        j = 1;
        num = 1;
        while j <= length(link(i).point)
            xx(num) = x(j);
            yy(num) = y(j);
            zz(num) = z(j);
            num = num + 1;
            j = j + 3;
        end
        if xx(num - 1) ~= x(length(link(i).point)) || yy(num - 1) ~= y(length(link(i).point)) || zz(num - 1) ~= z(length(link(i).point))
            xx(num) = x(length(link(i).point));
            yy(num) = y(length(link(i).point));
            zz(num) = z(length(link(i).point));
        end
        if length(xx) >= 2
            X = spline(linspace(0,1,length(xx)),[xx;yy;zz],linspace(0,1,20 * length(xx)));
            x_new = X(1,:);
            y_new = X(2,:);
            z_new = X(3,:);
        else
            x_new = xx;
            y_new = yy;
            z_new = zz;
        end
        %plot3(x_new,y_new,z_new,'linewidth',0.5,'color','k')
        curve = [x_new;y_new;z_new];%表示函数
        h = 0.01;%定义步长
        curve_1 = gradient(curve)./h;%求一阶导
        for ii = 1:length(x_new)
            if ii ~= length(x_new)
                sum_length = sum_length + ((x_new(ii+1)-x_new(ii))^2+(y_new(ii+1)-y_new(ii))^2+(z_new(ii+1)-z_new(ii))^2)^(1/2);
            end
            if round(x_new(ii)) ~= temp_x || round(y_new(ii)) ~= temp_y || round(z_new(ii)) ~= temp_z
                temp_x = round(x_new(ii));
                temp_y = round(y_new(ii));
                temp_z = round(z_new(ii));
            end
        end
    end
end

