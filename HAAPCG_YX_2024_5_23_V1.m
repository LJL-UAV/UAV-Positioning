clc; clear; warning('off'); close all;

baseFolder = 'F:\PC2-Data\UAV_GE1\Pairs_new4\YX\U2'; % 基础文件夹路径
folders = dir(fullfile(baseFolder, '*')); % 获取所有文件和文件夹
folders = folders([folders.isdir]); % 保留文件夹条目
folders = folders(~ismember({folders.name}, {'.', '..'})); % 移除'.'和'..'条目

patch_size = 128;

for i = 1:length(folders)
    folderName = folders(i).name; % 当前处理的文件夹名
    folderPath = fullfile(baseFolder, folderName); % 当前文件夹的完整路径

    % 指定图像文件路径
    str1 = fullfile(folderPath, '1360.png'); % image pair
    str2 = fullfile(folderPath, 'Ref.jpg'); % image pair

    if exist(str1, 'file') && exist(str2, 'file')
        image1 = im2uint8(imread(str1));
        image2 = im2uint8(imread(str2));
        image_1 = image1(:,:,1);
        image_2 = image2(:,:,1);

        % 2  初始参数设置
        K_weight=3;                        % 各向异性力矩图的加权值，处于（1~10），默认设置：3
        Max=3;                               % Number of levels in scale space，默认设置：3
        threshold = 0.3;                  % 特征点提取阈值，对SAR影像/强度图配色时，设置为：0.3；一般默认设置为：0.4
        scale_value=2;                  % 尺度缩放比例值，默认设置：1.6
        Path_Block=42;                   % 描述子邻域窗口大小， 默认设置：42；当需要更多特征点时，可以调大窗口。
        
        %%% 3 各向异性尺度空间
        [nonelinear_space_1]=HAPCG_nonelinear_space(image_1,Max,scale_value);
        [nonelinear_space_2]=HAPCG_nonelinear_space(image_2,Max,scale_value);

        %% 4  构建加权各向异性力矩图和相位一致性梯度及方向 
        [harris_function_1,gradient_1,angle_1]=HAPCG_Gradient_Feature(nonelinear_space_1,Max,K_weight);
        [harris_function_2,gradient_2,angle_2]=HAPCG_Gradient_Feature(nonelinear_space_2,Max,K_weight);

        % 5 特征检测,限制在5000
        position_1=Harris_extreme(harris_function_1,gradient_1,angle_1,Max,threshold);

        % 随机化 position_1 的行
        nRows = size(position_1, 1);  % 获取 position_1 的行数
        randIdx = randperm(nRows);    % 生成一个随机索引排列
        position_1 = position_1(randIdx, :);  % 重新排序 position_1
        % 检查是否有足够的行并提取前5000行
        if nRows >= 5000
            position_1 = position_1(1:5000, :);
        end

        position_2=Harris_extreme(harris_function_2,gradient_2,angle_2,Max,threshold);
        % 随机化 position_2 的行
        nRows = size(position_2, 1);  % 获取 position_2 的行数
        randIdx = randperm(nRows);    % 生成一个随机索引排列
        position_2 = position_2(randIdx, :);  % 重新排序 position_2
        
        % 检查是否有足够的行并提取前5000行
        if nRows >= 5000
            position_2 = position_2(1:5000, :);
        end

        % 写入kpts1到文件
        filenameKpts1 = fullfile(folderPath, 'HAPCG_kpts1.txt');
        dlmwrite(filenameKpts1, position_1(:, 1:2), 'delimiter', ' ', 'precision', '%f');
        % 写入kpts2到文件
        filenameKpts2 = fullfile(folderPath, 'HAPCG_kpts2.txt');
        dlmwrite(filenameKpts2, position_2(:, 1:2), 'delimiter', ' ', 'precision', '%f');


        % 6 特征描述
        descriptors_1=HAPCG_Logpolar_descriptors(gradient_1,angle_1,position_1,Path_Block);                                     
        descriptors_2=HAPCG_Logpolar_descriptors(gradient_2,angle_2,position_2,Path_Block); 

        % 将特征描述的点位置写入到文件
        filenameKpts1 = fullfile(folderPath, 'HAPCG_locs1.txt');
        dlmwrite(filenameKpts1, descriptors_1.locs(:, 1:4), 'delimiter', ' ', 'precision', '%f');
        filenameKpts2 = fullfile(folderPath, 'HAPCG_locs2.txt');
        dlmwrite(filenameKpts2, descriptors_2.locs(:,1:4), 'delimiter', ' ', 'precision', '%f');

        % 将特征描述的特征向量写入到文件
        filenameKpts1 = fullfile(folderPath, 'HAPCG_des1.txt');
        dlmwrite(filenameKpts1, descriptors_1.des(:,1:248), 'delimiter', ' ', 'precision', '%f');
        filenameKpts2 = fullfile(folderPath, 'HAPCG_des2.txt');
        dlmwrite(filenameKpts2, descriptors_2.des(:,1:248), 'delimiter', ' ', 'precision', '%f');


        % 特征匹配.
        [indexPairs,~] = matchFeatures(descriptors_1.des,descriptors_2.des,'MaxRatio',1,'MatchThreshold', 100);
        matchedPoints_1 = descriptors_1.locs(indexPairs(:, 1), :);
        matchedPoints_2 = descriptors_2.locs(indexPairs(:, 2), :);


        % 粗差剔除
        % 使用FSC函数计算两个特征点集之间的相似性变换矩阵H。
        [H,rmse]=FSC(matchedPoints_1,matchedPoints_2,'similarity',3);
        Y_=H*[matchedPoints_1(:,[1,2])';ones(1,size(matchedPoints_1,1))];
        Y_(1,:)=Y_(1,:)./Y_(3,:);
        Y_(2,:)=Y_(2,:)./Y_(3,:);
        E=sqrt(sum((Y_(1:2,:)-matchedPoints_2(:,[1,2])').^2));
        inliersIndex=E < 3;
        clearedPoints1 = matchedPoints_1(inliersIndex, :);
        clearedPoints2 = matchedPoints_2(inliersIndex, :);
        uni1=[clearedPoints1(:,[1,2]),clearedPoints2(:,[1,2])];
        [~,i,~]=unique(uni1,'rows','first');
        inliersPoints1=clearedPoints1(sort(i)',:);
        inliersPoints2=clearedPoints2(sort(i)',:);
        [inliersPoints_1,inliersPoints_2] = BackProjection(inliersPoints1,inliersPoints2,scale_value);

        %% 保存匹配图像
        folderPath = 'F:\PC2-Data\UAV_GE1\Pairs_new4\results\YX\HAPCG'; 
        % 检查文件夹是否存在，如果不存在则创建
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        % 创建文件名，包括其路径
        filename = fullfile(folderPath, sprintf('%s_HAPCG.png', folderName));
        % 保存当前图形到指定文件
        saveas(gcf, filename);

        %% 保存匹配点数量
        filename = fullfile(folderPath, 'HAPCG_NCM.txt');
        fid = fopen(filename, 'a');
        fprintf(fid, '%d,%d\n', size(clearedPoints1, 1), size(finalPoints1, 1));
        fclose(fid);
        close all;
    else
        fprintf('Files not found in %s. Skipping...\n', folderPath);
    end
end
