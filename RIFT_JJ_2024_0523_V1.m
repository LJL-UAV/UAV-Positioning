clc; clear; warning('off'); close all;

baseFolder = 'F:\PC2-Data\UAV_GE1\Pairs_new4\JJ\U2'; % 基础文件夹路径
folders = dir(fullfile(baseFolder, '*')); % 获取所有文件和文件夹
folders = folders([folders.isdir]); % 保留文件夹条目
folders = folders(~ismember({folders.name}, {'.', '..'})); % 移除'.'和'..'条目

for i = 1:length(folders)
    folderName = folders(i).name; % 当前处理的文件夹名
    folderPath = fullfile(baseFolder, folderName); % 当前文件夹的完整路径

    % 指定图像文件路径
    str1 = fullfile(folderPath, '1360.png'); % image pair
    str2 = fullfile(folderPath, 'Ref.jpg'); % image pair

    if exist(str1, 'file') && exist(str2, 'file')
        image1 = im2uint8(imread(str1));
        image2 = im2uint8(imread(str2));

        im1 = image1(:,:,1);
        im2 = image2(:,:,1);
        if size(im1,3) == 3
            im1 = rgb2gray(im1);
        end
        if size(im2,3) == 3
            im2 = rgb2gray(im2);
        end

        % 特征检测与描述
        [des_m1,des_m2] = RIFT_no_rotation_invariance(im1,im2,4,6,128);

        %% 特征点位置写到文件
        filenameKpts1 = fullfile(folderPath, 'RIFT_kpts1.txt');
        dlmwrite(filenameKpts1, des_m1.kps(:, 1:2), 'delimiter', ' ', 'precision', '%f');
        % 写入kpts2到文件
        filenameKpts2 = fullfile(folderPath, 'RIFT_kpts2.txt');
        dlmwrite(filenameKpts2, des_m2.kps(:, 1:2), 'delimiter', ' ', 'precision', '%f');
        
        %% 特征点描述写到文件
        filenameKpts1 = fullfile(folderPath, 'RIFT_des1.txt');
        dlmwrite(filenameKpts1, des_m1.des(:, 1:216), 'delimiter', ' ', 'precision', '%f');
        % 写入kpts2到文件
        filenameKpts2 = fullfile(folderPath, 'RIFT_des2.txt');
        dlmwrite(filenameKpts2, des_m2.des(:, 1:216), 'delimiter', ' ', 'precision', '%f');


        % 特征匹配.
        [indexPairs,matchmetric] = matchFeatures(des_m1.des,des_m2.des,'MaxRatio',1,'MatchThreshold', 100);
        matchedPoints1 = des_m1.kps(indexPairs(:, 1), :);
        matchedPoints2 = des_m2.kps(indexPairs(:, 2), :);
        [matchedPoints2,IA]=unique(matchedPoints2,'rows');
        matchedPoints1=matchedPoints1(IA,:);
        
        % 根据匹配对的索引，从特征点矩阵中提取匹配点的位置信息。
        matchedPoints1 = des_m1.kps(indexPairs(:, 1), 1:2);
        matchedPoints2 = des_m2.kps(indexPairs(:, 2), 1:2);
        
        % 它使用唯一性去除重复的匹配点。
        [matchedPoints2,IA]=unique(matchedPoints2,'rows');
        matchedPoints1=matchedPoints1(IA,:);
        
        % 粗差剔除
        % 使用FSC函数计算两个特征点集之间的相似性变换矩阵H。
        % matchedPoints1和matchedPoints2是匹配点的位置信息。'similarity'参数指定了相似性变换模型，3是指定的RANSAC算法的迭代次数。
        H=FSC(matchedPoints1,matchedPoints2,'similarity',3);
        
        % 应用相似性变换矩阵H将matchedPoints1中的点变换到matchedPoints2的坐标系中。变换后的结果存储在矩阵Y_中，然后对Y_进行归一化处理。
        Y_=H*[matchedPoints1';ones(1,size(matchedPoints1,1))];
        Y_(1,:)=Y_(1,:)./Y_(3,:);
        Y_(2,:)=Y_(2,:)./Y_(3,:);
        
        %计算变换后的点与原始matchedPoints2的欧氏距离，并根据阈值3将距离小于阈值的点标记为内点。
        % 然后，根据内点的索引，将内点存储在cleanedPoints1和cleanedPoints2中。
        E=sqrt(sum((Y_(1:2,:)-matchedPoints2').^2));
        inliersIndex=E<3;
        cleanedPoints1 = matchedPoints1(inliersIndex, :);
        cleanedPoints2 = matchedPoints2(inliersIndex, :);
        [cleanedPoints2,IA] = unique(cleanedPoints2,'rows');
        cleanedPoints1 = cleanedPoints1(IA,:);
       

        %% 保存匹配图像
        folderPath = 'F:\PC2-Data\UAV_GE1\Pairs_new4\results\JJ\RIFT';  % 你可以将其修改为任何具体的文件夹路径，例如 'C:\Users\YourUsername\Documents\Matlab\a'
        % 检查文件夹是否存在，如果不存在则创建
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        % 创建文件名，包括其路径
        filename = fullfile(folderPath, sprintf('%s_RIFT.png', folderName));
        % 保存当前图形到指定文件
        saveas(gcf, filename);

        %% 保存匹配点数量
        filename = fullfile(folderPath, 'RIFT_NCM.txt');
        fid = fopen(filename, 'a');
        fprintf(fid, '%d,%d\n', size(cleanedPoints1, 1), size(finalPoints1, 1));
        fclose(fid);
        close all;
        
    else
        fprintf('Files not found in %s. Skipping...\n', folderPath);
    end
end
