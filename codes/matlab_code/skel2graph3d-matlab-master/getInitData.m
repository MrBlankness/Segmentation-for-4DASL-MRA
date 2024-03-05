clear;
clc;

load A003_Result_maskB_Filter.mat
mask = img;
skel = Skeleton3D(imbinarize(mask));
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
[~,node,link] = Skel2Graph3D(skel,0);
rData = img * 0;

for i = 1:w
    for j = 1:l
        for k = 1:h
            if(skel(i,j,k)~=0)
                x = i;
                y = j;
                z = k;
                r = 1;
                while(1)
                    v = 0;
                    vv = 0;
                    for xx = x-r:x+r
                        if(xx<1||xx>w)
                            continue;
                        end
                        for yy = y-r:y+r
                            if(yy<1||yy>l)
                                continue;
                            end
                            for zz = z-r:z+r
                                if(zz<1||zz>h)
                                    continue;
                                end
                                if(sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz))<=r)
                                    v = v + 1;
                                    if(mask(xx,yy,zz)==1)
                                        vv = vv +1;
                                    end
                                end
                            end
                        end
                    end
                    if(v>vv)
                        % fprintf('x=%d,y=%d,z=%d,v=%d,vv=%d,r=%d\n',x,y,z,v,vv,r-1);
                        rData(x,y,z) = r;
                        break;
                    else
                        r = r + 1;
                    end
                end
            end
        end
    end
end


testData = img * 0;
for i = 1:w
    for j = 1:l
        for k = 1:h
            if(rData(i,j,k)~=0)
                x = i;
                y = j;
                z = k;
                r = rData(x,y,z);
                for xx = x-r:x+r
                    if(xx<1||xx>w)
                        continue;
                    end
                    for yy = y-r:y+r
                        if(yy<1||yy>l)
                            continue;
                        end
                        for zz = z-r:z+r
                            if(zz<1||zz>h)
                                continue;
                            end
                            if(sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz))<=r)
                                testData(xx,yy,zz) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

save('skel_003.mat','skel');
save('testData_003.mat','testData');