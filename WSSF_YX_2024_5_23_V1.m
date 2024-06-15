%%
% -----！！！！-----其分别监测角点与斑点，然后分别进行特征描述与匹配，最后再将匹配结果进行合并。----！！！！---
%%
clc; 
clear; 
warning('off');
close all;
beep off;
addpath(genpath('PSATF'));
addpath(genpath('Others'));


baseFolder = 'F:\PC2-Data\UAV_GE1\Pairs_new4\YX\U2'; % 基础文件夹路径
folders = dir(fullfile(baseFolder, '*')); % 获取所有文件和文件夹
folders = folders([folders.isdir]); % 保留文件夹条目
folders = folders(~ismember({folders.name}, {'.', '..'})); % 移除'.'和'..'条目

%% 2  参数设置
Path_Block=48;                   
sigma_1=1.6;   
ratio=2^(1/3);                     
ScaleValue = 1.6;
nOctaves = 3;
filter = 5;
Scale ='YES';

for i = 1:length(folders)
    folderName = folders(i).name; % 当前处理的文件夹名
    folderPath = fullfile(baseFolder, folderName); % 当前文件夹的完整路径

    % 指定图像文件路径
    str1 = fullfile(folderPath, '1360.png'); % image pair
    str2 = fullfile(folderPath, 'Ref.jpg'); % image pair

    if exist(str1, 'file') && exist(str2, 'file')
        %
        image_3 = im2uint8(imread(str1));
        image_4 = im2uint8(imread(str2));
        image_1 = im2double(image_3);
        image_2 = im2double(image_4);

        % 构造影像空间
        t1=clock;
        tic;
        [nonelinear_space_1,E_space_1,Max_space_1,Min_space_1,Phase_space_1]=Create_Image_space(image_1,nOctaves,Scale, ScaleValue, ratio,sigma_1,filter);
        [nonelinear_space_2,E_space_2,Max_space_2,Min_space_2,Phase_space_2]=Create_Image_space(image_2,nOctaves,Scale, ScaleValue, ratio,sigma_1,filter);
        disp(['构造影像尺度空间花费时间：',num2str(toc),' 秒']);

        % 特征检测
        tic;
        [Bolb_KeyPts_1,Corner_KeyPts_1,Bolb_gradient_1,Corner_gradient_1,Bolb_angle_1,Corner_angle_1]  =  WSSF_features(nonelinear_space_1,E_space_1,Max_space_1,Min_space_1,Phase_space_1,sigma_1,ratio,Scale,nOctaves);
        [Bolb_KeyPts_2,Corner_KeyPts_2,Bolb_gradient_2,Corner_gradient_2,Bolb_angle_2,Corner_angle_2]  =  WSSF_features(nonelinear_space_2,E_space_2,Max_space_2,Min_space_2,Phase_space_2,sigma_1,ratio,Scale,nOctaves);
        disp(['特征点提取花费时间:  ',num2str(toc),' 秒']);

        %% 特征点位置写到文件
        % 写入kpts1到文件
        filenameKpts1 = fullfile(folderPath, 'WSSF_kpts1.txt');
        % 以覆盖模式写入第一组数据，保证文件是新的或清空原有内容
        dlmwrite(filenameKpts1, Bolb_KeyPts_1(:,1:2), 'delimiter', ' ', 'precision', '%f');
        % 以追加模式写入第二组数据
        dlmwrite(filenameKpts1, Corner_KeyPts_1(:,1:2), 'delimiter', ' ', 'precision', '%f', '-append');
        
        % 写入kpts2到文件
        filenameKpts2 = fullfile(folderPath, 'WSSF_kpts2.txt');
        % 以覆盖模式写入第一组数据
        dlmwrite(filenameKpts2, Bolb_KeyPts_2(:,1:2), 'delimiter', ' ', 'precision', '%f');
        % 以追加模式写入第二组数据
        dlmwrite(filenameKpts2, Corner_KeyPts_2(:,1:2), 'delimiter', ' ', 'precision', '%f', '-append');


        % 特征描述
        tic;
        Bolb_descriptors_1 = GLOH_descriptors(Bolb_gradient_1, Bolb_angle_1, Bolb_KeyPts_1, Path_Block, ratio,sigma_1);
        Corner_descriptors_1 = GLOH_descriptors(Corner_gradient_1, Corner_angle_1, Corner_KeyPts_1, Path_Block, ratio,sigma_1);
        Bolb_descriptors_2 = GLOH_descriptors(Bolb_gradient_2, Bolb_angle_2, Bolb_KeyPts_2, Path_Block, ratio,sigma_1);
        Corner_descriptors_2 = GLOH_descriptors(Corner_gradient_2, Corner_angle_2, Corner_KeyPts_2, Path_Block, ratio,sigma_1);
        disp(['特征描述子花费时间:  ',num2str(toc),' 秒']); 
            
        %% 特征描述写到文件
        filenameKpts1 = fullfile(folderPath, 'WSSF_des1.txt');
        % 以覆盖模式写入第一组数据，保证文件是新的或清空原有内容
        dlmwrite(filenameKpts1, Bolb_descriptors_1.des(:,1:396), 'delimiter', ' ', 'precision', '%f');
        % 以追加模式写入第二组数据
        dlmwrite(filenameKpts1, Corner_descriptors_1.des(:,1:396), 'delimiter', ' ', 'precision', '%f', '-append');
        
        % 写入kpts2到文件
        filenameKpts2 = fullfile(folderPath, 'WSSF_des2.txt');
        % 以覆盖模式写入第一组数据
        dlmwrite(filenameKpts2, Bolb_descriptors_2.des(:,1:396), 'delimiter', ' ', 'precision', '%f');
        % 以追加模式写入第二组数据
        dlmwrite(filenameKpts2, Corner_descriptors_2.des(:,1:396), 'delimiter', ' ', 'precision', '%f', '-append');


        % 特征匹配与交叉粗提纯
        [indexPairs,~]= matchFeatures(Bolb_descriptors_1.des,Bolb_descriptors_2.des,'MaxRatio',1,'MatchThreshold', 50,'Unique',true ); 
        [matchedPoints_1_1,matchedPoints_1_2] = BackProjection(Bolb_descriptors_1.locs(indexPairs(:, 1), :),Bolb_descriptors_2.locs(indexPairs(:, 2), :),ScaleValue); 
        [indexPairs,~]= matchFeatures(Corner_descriptors_1.des,Corner_descriptors_2.des,'MaxRatio',1,'MatchThreshold', 50,'Unique',true ); 
        [matchedPoints_2_1,matchedPoints_2_2] = BackProjection(Corner_descriptors_1.locs(indexPairs(:, 1), :),Corner_descriptors_2.locs(indexPairs(:, 2), :),ScaleValue); 
        
        matchedPoints_1 = [matchedPoints_1_1;matchedPoints_2_1];
        matchedPoints_2 = [matchedPoints_1_2;matchedPoints_2_2];
        allNCM = size(matchedPoints_1,1);
        
        [H1,rmse]=FSC(matchedPoints_1,matchedPoints_2,'similarity',3);
        Y_=H1*[matchedPoints_1(:,[1,2])';ones(1,size(matchedPoints_1,1))];
        Y_(1,:)=Y_(1,:)./Y_(3,:);
        Y_(2,:)=Y_(2,:)./Y_(3,:);
        E=sqrt(sum((Y_(1:2,:)-matchedPoints_2(:,[1,2])').^2));
        inliersIndex=E < 3;
        cleanedPoints1 = matchedPoints_1(inliersIndex, :);
        cleanedPoints2 = matchedPoints_2(inliersIndex, :);
        [cleanedPoints2,IA]=unique(cleanedPoints2,'rows');
        cleanedPoints1=cleanedPoints1(IA,:);

        %% 保存匹配图像
        folderPath = 'F:\PC2-Data\UAV_GE1\Pairs_new4\results\YX\WSSF'; 
        % 检查文件夹是否存在，如果不存在则创建
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        % 创建文件名，包括其路径
        filename = fullfile(folderPath, sprintf('%s_WSSF.png', folderName));
        % 保存当前图形到指定文件
        saveas(gcf, filename);

        %% 保存匹配点数量
        filename = fullfile(folderPath, 'WSSF_NCM.txt');
        fid = fopen(filename, 'a');
        fprintf(fid, '%d,%d\n', size(cleanedPoints1, 1), size(finalPoints1, 1));
        fclose(fid);
        close all;
    else
        fprintf('Files not found in %s. Skipping...\n', folderPath);
    end
end
