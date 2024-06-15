clc; clear; warning('off'); close all;

baseFolder = 'F:\PC2-Data\UAV_GE1\Pairs_new4\JJ\U2'; % �����ļ���·��
folders = dir(fullfile(baseFolder, '*')); % ��ȡ�����ļ����ļ���
folders = folders([folders.isdir]); % �����ļ�����Ŀ
folders = folders(~ismember({folders.name}, {'.', '..'})); % �Ƴ�'.'��'..'��Ŀ

for i = 1:length(folders)
    folderName = folders(i).name; % ��ǰ������ļ�����
    folderPath = fullfile(baseFolder, folderName); % ��ǰ�ļ��е�����·��

    % ָ��ͼ���ļ�·��
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

        % �������������
        [des_m1,des_m2] = RIFT_no_rotation_invariance(im1,im2,4,6,128);

        %% ������λ��д���ļ�
        filenameKpts1 = fullfile(folderPath, 'RIFT_kpts1.txt');
        dlmwrite(filenameKpts1, des_m1.kps(:, 1:2), 'delimiter', ' ', 'precision', '%f');
        % д��kpts2���ļ�
        filenameKpts2 = fullfile(folderPath, 'RIFT_kpts2.txt');
        dlmwrite(filenameKpts2, des_m2.kps(:, 1:2), 'delimiter', ' ', 'precision', '%f');
        
        %% ����������д���ļ�
        filenameKpts1 = fullfile(folderPath, 'RIFT_des1.txt');
        dlmwrite(filenameKpts1, des_m1.des(:, 1:216), 'delimiter', ' ', 'precision', '%f');
        % д��kpts2���ļ�
        filenameKpts2 = fullfile(folderPath, 'RIFT_des2.txt');
        dlmwrite(filenameKpts2, des_m2.des(:, 1:216), 'delimiter', ' ', 'precision', '%f');


        % ����ƥ��.
        [indexPairs,matchmetric] = matchFeatures(des_m1.des,des_m2.des,'MaxRatio',1,'MatchThreshold', 100);
        matchedPoints1 = des_m1.kps(indexPairs(:, 1), :);
        matchedPoints2 = des_m2.kps(indexPairs(:, 2), :);
        [matchedPoints2,IA]=unique(matchedPoints2,'rows');
        matchedPoints1=matchedPoints1(IA,:);
        
        % ����ƥ��Ե����������������������ȡƥ����λ����Ϣ��
        matchedPoints1 = des_m1.kps(indexPairs(:, 1), 1:2);
        matchedPoints2 = des_m2.kps(indexPairs(:, 2), 1:2);
        
        % ��ʹ��Ψһ��ȥ���ظ���ƥ��㡣
        [matchedPoints2,IA]=unique(matchedPoints2,'rows');
        matchedPoints1=matchedPoints1(IA,:);
        
        % �ֲ��޳�
        % ʹ��FSC�����������������㼯֮��������Ա任����H��
        % matchedPoints1��matchedPoints2��ƥ����λ����Ϣ��'similarity'����ָ���������Ա任ģ�ͣ�3��ָ����RANSAC�㷨�ĵ���������
        H=FSC(matchedPoints1,matchedPoints2,'similarity',3);
        
        % Ӧ�������Ա任����H��matchedPoints1�еĵ�任��matchedPoints2������ϵ�С��任��Ľ���洢�ھ���Y_�У�Ȼ���Y_���й�һ������
        Y_=H*[matchedPoints1';ones(1,size(matchedPoints1,1))];
        Y_(1,:)=Y_(1,:)./Y_(3,:);
        Y_(2,:)=Y_(2,:)./Y_(3,:);
        
        %����任��ĵ���ԭʼmatchedPoints2��ŷ�Ͼ��룬��������ֵ3������С����ֵ�ĵ���Ϊ�ڵ㡣
        % Ȼ�󣬸����ڵ�����������ڵ�洢��cleanedPoints1��cleanedPoints2�С�
        E=sqrt(sum((Y_(1:2,:)-matchedPoints2').^2));
        inliersIndex=E<3;
        cleanedPoints1 = matchedPoints1(inliersIndex, :);
        cleanedPoints2 = matchedPoints2(inliersIndex, :);
        [cleanedPoints2,IA] = unique(cleanedPoints2,'rows');
        cleanedPoints1 = cleanedPoints1(IA,:);
       

        %% ����ƥ��ͼ��
        folderPath = 'F:\PC2-Data\UAV_GE1\Pairs_new4\results\JJ\RIFT';  % ����Խ����޸�Ϊ�κξ�����ļ���·�������� 'C:\Users\YourUsername\Documents\Matlab\a'
        % ����ļ����Ƿ���ڣ�����������򴴽�
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        % �����ļ�����������·��
        filename = fullfile(folderPath, sprintf('%s_RIFT.png', folderName));
        % ���浱ǰͼ�ε�ָ���ļ�
        saveas(gcf, filename);

        %% ����ƥ�������
        filename = fullfile(folderPath, 'RIFT_NCM.txt');
        fid = fopen(filename, 'a');
        fprintf(fid, '%d,%d\n', size(cleanedPoints1, 1), size(finalPoints1, 1));
        fclose(fid);
        close all;
        
    else
        fprintf('Files not found in %s. Skipping...\n', folderPath);
    end
end
