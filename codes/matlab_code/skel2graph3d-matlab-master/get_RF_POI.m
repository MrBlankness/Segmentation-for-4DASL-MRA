clear;
clc;

load A003_Result_maskB_Filter.mat
mask = img;
skel = Skeleton3D(imbinarize(mask));
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
[~,node,link] = Skel2Graph3D(skel,0);
rData = mask * 0;

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

poiData = mask * 0;
for i=1:length(node)
    for j=1:length(node(i).links)
        for k=2:length(link(node(i).links(j)).point)-1
            [x0,y0,z0]=ind2sub([w,l,h],link(node(i).links(j)).point(k-1));
            [x,y,z]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
            [x1,y1,z1]=ind2sub([w,l,h],link(node(i).links(j)).point(k+1));
            dX = [-1,0,1];
            dY = [rData(x0,y0,z0),rData(x,y,z),rData(x1,y1,z1)];
            A = polyfit(dX,dY,2);
            a = A(1);
            b = A(2);
            c = A(3);
            if(1/(2*abs(a))>1 && 1/(2*abs(a))<4*dY(2))
                r = rData(x,y,z);
                for ii = x-r:x+r
                    if(ii<1||ii>w)
                        continue
                    end
                    for jj = y-r:y+r
                        if(jj<1||jj>l)
                            continue
                        end
                        for kk = z-r:z+r
                            if(kk<1||kk>h)
                                continue
                            else
                                if(sqrt((x-ii)^2+(y-jj)^2+(z-kk)^2)<=r&&sqrt((x-ii)^2+(y-jj)^2+(z-kk)^2)>(r-1))
                                    poiData(ii,jj,kk) = 1;
                                end
                            end
                        end
                    end
                end
%                 judgeData(x,y,z) = 1;
%                 number2 = number2 + 1;
            end
        end
    end
end

num = 0;
for i=1:length(node)
    for j=1:length(node(i).links)
        testData = mask * 0;
        maxR = 0;
        for k=1:length(link(node(i).links(j)).point)
            [x,y,z]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
            r = rData(x,y,z);
            if(r > maxR)
                maxR = r;
            end
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
        poiPartData = poiData .* testData;
        for ii = 1:w
            for jj = 1:l
                for kk = 1:h
                    if(poiPartData(ii,jj,kk)~=0)
                        num = num + 1;
                        poiResData(num).x = ii;
                        poiResData(num).y = jj;
                        poiResData(num).z = kk;
                        poiResData(num).d = 10000;
                        if(length(link(node(i).links(j)).point) > 4 * maxR)
                            poiResData(num).isLongTrunk = 1;
                        else
                            poiResData(num).isLongTrunk = 0;
                        end
                        for k=1:length(link(node(i).links(j)).point)
                            [x,y,z]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
                            dis = sqrt((x-ii)*(x-ii)+(y-jj)*(y-jj)+(z-kk)*(z-kk));
                            if(dis<poiResData(num).d)
                                poiResData(num).centerX = x;
                                poiResData(num).centerY = y;
                                poiResData(num).centerZ = z;
                                poiResData(num).d = dis;
                            end
                        end
                        if(poiResData(num).d == 10000)
                            poiResData(num).d = 0;
                        end
                        poiResData(num).r = rData(poiResData(num).centerX,poiResData(num).centerY,poiResData(num).centerZ);
                        poiResData(num).Planeness = 1;
                        for mm = ii-1:ii+1
                            if(mm > w || mm < 1)
                                continue;
                            end
                            for nn = jj-1:jj+1
                                if(nn > l || nn < 1)
                                    continue;
                                end
                                for pp = kk-1:kk+1
                                    if( pp > h || pp < 1)
                                        continue;
                                    end
                                    if(testData(mm,nn,pp)~=0)
                                        xx0 = mm;
                                        yy0 = nn;
                                        zz0 = pp;
                                        for mmm = ii-1:ii+1
                                            if(mmm > w || mmm < 1)
                                                continue;
                                            end
                                            for nnn = jj-1:jj+1
                                                if(nnn > l || nnn < 1)
                                                    continue;
                                                end
                                                for ppp = kk-1:kk+1
                                                    if( ppp > h || ppp < 1)
                                                        continue;
                                                    end
                                                    if(testData(mmm,nnn,ppp)~=0)
                                                        xx1 = mmm;
                                                        yy1 = nnn;
                                                        zz1 = ppp;
                                                        if(xx0==xx1&&yy0==yy1&&zz0==zz1)
                                                            continue;
                                                        else
                                                            A = [xx0-ii,yy0-jj,zz0-kk];
                                                            B = [xx1-ii,yy1-jj,zz1-kk];
                                                            Planeness = dot(A,B)/norm(A)/norm(B);
                                                            if(Planeness < poiResData(num).Planeness)
                                                                poiResData(num).Planeness = Planeness;
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        rr = poiResData(num).r;
                        v = 0;
                        for mm = ii - rr : ii + rr
                            if(mm > w || mm < 1)
                                continue;
                            end
                            for nn = jj - rr : jj + rr
                                if(nn > l || nn < 1)
                                    continue;
                                end
                                for pp = kk - rr : kk + rr
                                    if(pp > h || pp < 1)
                                        continue;
                                    end
                                    if(sqrt((mm-ii)*(mm-ii)+(nn-jj)*(nn-jj)+(pp-kk)*(pp-kk))<rr&&testData(mm,nn,pp)~=0)
                                        v = v + 1;
                                    end
                                end
                            end
                        end
                        poiResData(num).CS = v / (4/3*pi*(rr^3));
                    end
                end
            end
        end
    end
end


num = 0;
for i = 1:length(poiResData)
    if(poiResData(i).r==1)
        num = num + 1;
        index(num) = i;
        continue
    end
    if(poiResData(i).centerX==poiResData(i).x&&poiResData(i).centerY==poiResData(i).y&&poiResData(i).centerZ==poiResData(i).z)
        num = num + 1;
        index(num) = i;
    end
end

poiResData(index) = [];
save('poiResData_RF_003.mat','poiResData');