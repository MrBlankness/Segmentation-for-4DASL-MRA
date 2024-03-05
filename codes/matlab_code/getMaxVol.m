function data = getMaxVol(data)
    D = bwconncomp(data);
    max_voxel = 0;
    index = 0;
    for i = 1:length(D.PixelIdxList)
        if length(D.PixelIdxList{1,i}) > max_voxel
            max_voxel = length(D.PixelIdxList{1,i});
            index = i;
        end
    end
    D = bwconncomp(data);
    for i = 1:length(D.PixelIdxList)
        if i ~= index
            data(D.PixelIdxList{1,i}) = 0;
        end
    end
end

