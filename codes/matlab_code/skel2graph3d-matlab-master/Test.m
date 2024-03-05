clear;
clc;
% fprintf('start\n');
% load A133_Result_maskB.mat

% figure(),volshow(img);
% disData = img * 0;
% w = size(disData,1);
% l = size(disData,2);
% h = size(disData,3);
% i=128;
% j=128;
% k=0;
% for ii = i-5:i+5
%     if(ii<1||ii>w)
%         continue;
%     end
%     for jj = j-5:j+5
%         if(jj<1||jj>l)
%             continue;
%         end
%         for kk = k-5:k+5
%             if(kk<1||k>h)
%                 continue;
%             end
%             disData(ii,jj,kk)=1;
%         end
%     end
% end
%             
% 
% figure(),volshow(disData);
% 
% figure('Name','原图'),volshow(img);



load A133_Result_maskB_Filter.mat
mask = img;
skel = Skeleton3D(imbinarize(mask));
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
[~,node,link] = Skel2Graph3D(skel,0);

rData = img*0;

% 提取半径
for i=1:length(node)    
    for j=1:length(node(i).links)
        for k=1:length(link(node(i).links(j)).point)
            [x,y,z]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
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

% 根据半径重建三维数据
% testData = img * 0;
% for i=1:length(node)
%     for j=1:length(node(i).links)
%         for k=1:length(link(node(i).links(j)).point)
%             [x,y,z]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
%             r = rData(x,y,z);
%             for xx = x-r:x+r
%                 if(xx<1||xx>w)
%                     continue;
%                 end
%                 for yy = y-r:y+r
%                     if(yy<1||yy>l)
%                         continue;
%                     end
%                     for zz = z-r:z+r
%                         if(zz<1||zz>h)
%                             continue;
%                         end
%                         if(sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz))<=r)
%                             testData(xx,yy,zz) = 1;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

% figure('Name','过滤'),volshow(img);
% figure('Name','重建后'),volshow(testData);
% figure('Name','骨架'),volshow(skel);

% 去找最低的点作为距离变换的初始点
disData = img * 0;
judge = 0;
for i = 1:w
    for j = 1:l
        for k = 1:h
            if(skel(i,j,k) == 1)
                disData(i,j,k) = 1;
                judge = 1;
                break
            end
        end
        if(judge == 1)
            break;
        end
    end
    if(judge == 1)
        break;
    end
end
% volshow(disData);
disData = bwdist(disData);
disData = disData .* img;
resData = img * 0;
for i = 1:w
    for j = 1:l
        for k = 1:h
            if(disData(i,j,k)~=0)
                judge = 0;
                for ii = i-1:2:i+1
                    if(ii<1||ii>w)
                        continue;
                    end
                    for jj = j-1:2:j+1
                        if(jj<1||jj>l)
                            continue;
                        end
                        for kk = k-1:2:k+1
                            if(kk<1||kk>h)
                                continue;
                            end
                            if(disData(ii,jj,kk)>disData(i,j,k))
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
                if(judge==0)
                    resData(i,j,k)=1;
                end
            end
        end
    end
end
                    
figure(),volshow(resData);
save('resData.mat','resData');

maskPath = 'D:\idmWorkspace\Data Set 1 Digital subtraction angiography images and aneurysm mask and geometry data for tasks 1, 2, and 3_2\CADA-Training_MaskImages-NIFTI\A133_masks.nii.gz';
nii=load_untouch_nii(maskPath);
maskData=nii.img;
figure(),volshow(maskData);

% for i = 1:w
%     for j = 1:l
%         for k = 1:h
%             if(resData(i,j,k)==1&&maskData(i,j,k)==1)
%                 fprintf('x=%d,y=%d,z=%d\n',i,j,k);
%             end
%         end
%     end
% end
       


















% for i=1:length(node)
%     for j=1:length(node(i).links)
%         for k=2:length(link(node(i).links(j)).point)-1
%             [x0,y0,z0]=ind2sub([w,l,h],link(node(i).links(j)).point(k-1));
%             [x,y,z]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
%             [x1,y1,z1]=ind2sub([w,l,h],link(node(i).links(j)).point(k+1));
%             if(rData(x,y,z)>rData(x0,y0,z0) && rData(x,y,z)>rData(x1,y1,z1))
%                 fprintf('x=%d,y=%d,z=%d\n',x,y,z);
%             end
%         end
%     end
% end
% number1 = 0;
% judgeData = img*0;
% for i=1:length(node)
%     for j=1:length(node(i).links)
%         rrData = img * 0;
%         disData = img * 0;
%         for k=1:length(link(node(i).links(j)).point)
%             [x,y,z]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
%             if(k==1)
%                 disData(x,y,z) = 1;
%                 fprintf('bwdist start\n');
%                 disData = bwdist(disData);
%                 fprintf('bwdist end\n');
%             end
%             r = rData(x,y,z);
%             for xx = x-r:x+r
%                 if(xx<1||xx>w)
%                     continue
%                 end
%                 for yy = y-r:y+r
%                     if(yy<1||yy>l)
%                         continue
%                     end
%                     for zz = z-r:z+r
%                         if(zz<1||zz>h)
%                             continue
%                         end
%                         if(sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz))<=r)
%                             rrData(xx,yy,zz) = 1;
%                         end
%                     end
%                 end
%             end   
%         end
%         resData = disData.*rrData;
%         for ii = 1:w
%             for jj = 1:l
%                 for kk = 1:h
%                     judge = 0;
%                     for iii = ii - 1:ii + 1
%                         if(iii>w||iii<1)
%                             break;
%                         end
%                         for jjj = jj - 1:jj + 1
%                             if(jjj>l||jjj<1)
%                                 break;
%                             end
%                             for kkk = kk - 1:kk + 1
%                                 if(kkk>h||kkk<1)
%                                     break;
%                                 end
%                                 if(resData(ii,jj,kk)<resData(iii,jjj,kkk))
%                                     judge = 1;
%                                 end
%                             end
%                         end
%                     end
%                     if(judge==0)
%                         judgeData(ii,jj,kk)=1;
%                         number1 = number1 + 1;
%                     end
%                 end
%             end
%         end                     
%     end
% end
% fprintf('DT：%d\n',number1);
% maskPath = 'D:\idmWorkspace\Data Set 1 Digital subtraction angiography images and aneurysm mask and geometry data for tasks 1, 2, and 3_2\CADA-Training_MaskImages-NIFTI\A133_masks.nii.gz';
% nii=load_untouch_nii(maskPath);
% maskData=nii.img;
% for i = 1:w
%     for j = 1:l
%         for k = 1:h
%             if(judgeData(i,j,k)==1&&maskData(i,j,k)==1)
%                 fprintf('x=%d,y=%d,z=%d\n',i,j,k);
%             end
%         end
%     end
% end
% 
number2 = 0;
judgeData = img * 0;
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
            if(1/(2*abs(a))>1&&1/(2*abs(a))<4*dY(2))
                judgeData(x,y,z) = 1;
                number2 = number2 + 1;
            end
        end
    end
end
figure(),volshow(judgeData);
fprintf('RF：%d\n',number2);
for i = 1:w
    for j = 1:l
        for k = 1:h
            if(judgeData(i,j,k)==1&&maskData(i,j,k)==1)
                fprintf('x=%d,y=%d,z=%d\n',i,j,k);
            end
        end
    end
end