clear;
clc;

load A133_Result_maskB_Filter.mat
mask = img;
spacing = 0.542946994304657;
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
                        rData(x,y,z) = r - 1;
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

save('testData.mat','testData');
save('skel.mat','skel');


differenceData = mask - testData;
index = (differenceData == 1);
differenceData = img * 0;
differenceData(index) = 1;
D = bwconncomp(differenceData);


indexNum = 0;
for i = 1:length(D.PixelIdxList)
    if(length(D.PixelIdxList{1,i}) > 4/3*pi*(2/spacing))
        x = [];
        y = [];
        z = [];
        num = 0;
        temp = D.PixelIdxList{1,i};
        for j = 1:length(temp)
            num = num + 1;
            subarray = myind2sub(size(differenceData), temp(j,1));
            z(num) = subarray(3);
            y(num) = subarray(2);
            x(num) = subarray(1);
        end
        centerX = round(sum(x)/length(x));
        centerY = round(sum(y)/length(y));
        centerZ = round(sum(z)/length(z));
%         oneData = img * 0;
%         oneData(D.PixelIdxList{1,i}) = 1;
%         imgrayx = zeros(l,h);
%         imgrayy = zeros(w,h);
%         imgrayz = zeros(w,l);
%         for ii = 1:w
%             imgrayx(ii,:,:)=ii.*oneData(ii,:,:);
%         end
%         for jj = 1:l
%             imgrayy(:,jj,:)=jj.*oneData(:,jj,:);
%         end
%         for kk = 1:h
%             imgrayz(:,:,kk)=kk.*oneData(:,:,kk);
%         end
%         m = length(D.PixelIdxList{1,i});
%         centerX = sum(sum(sum(imgrayx)))/m;
%         centerY = sum(sum(sum(imgrayy)))/m;
%         centerZ = sum(sum(sum(imgrayz)))/m;
        rr = 1;
        judge = 0;
        while(1)
            for ii = centerX-r:centerX+r
                if(ii<1 || ii>w)
                    continue;
                end
                for jj = centerY-r:centerY+r
                    if(jj<1 || jj>l)
                        continue;
                    end
                    for kk = centerZ-r:centerZ+r
                        if(kk<1 || kk>h)
                            continue;
                        end
                        if(differenceData(ii,jj,kk)==0)
                            judge = 1;
                            break;
                        end
                    end
                    if(judge==1)
                        break;
                    end
                end
                if(judge==1)
                    break;
                end
            end
            if(judge==1)
                break;
            end
            rr = rr + 1;
        end
        rr = rr - 1;
        if(length(temp) > 9/8*r^3)
            indexNum = indexNum + 1;
            differencePOI(indexNum).x = centerX;
            differencePOI(indexNum).y = centerY;
            differencePOI(indexNum).z = centerZ;
        end
    end
end

