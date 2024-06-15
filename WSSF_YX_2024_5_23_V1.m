%%
% -----��������-----��ֱ���ǵ���ߵ㣬Ȼ��ֱ��������������ƥ�䣬����ٽ�ƥ�������кϲ���----��������---
%%
clc; 
clear; 
warning('off');
close all;
beep off;
addpath(genpath('PSATF'));
addpath(genpath('Others'));


baseFolder = 'F:\PC2-Data\UAV_GE1\Pairs_new4\YX\U2'; % �����ļ���·��
folders = dir(fullfile(baseFolder, '*')); % ��ȡ�����ļ����ļ���
folders = folders([folders.isdir]); % �����ļ�����Ŀ
folders = folders(~ismember({folders.name}, {'.', '..'})); % �Ƴ�'.'��'..'��Ŀ

%% 2  ��������
Path_Block=48;                   
sigma_1=1.6;   
ratio=2^(1/3);                     
ScaleValue = 1.6;
nOctaves = 3;
filter = 5;
Scale ='YES';

for i = 1:length(folders)
    folderName = folders(i).name; % ��ǰ������ļ�����
    folderPath = fullfile(baseFolder, folderName); % ��ǰ�ļ��е�����·��

    % ָ��ͼ���ļ�·��
    str1 = fullfile(folderPath, '1360.png'); % image pair
    str2 = fullfile(folderPath, 'Ref.jpg'); % image pair

    if exist(str1, 'file') && exist(str2, 'file')
        %
        image_3 = im2uint8(imread(str1));
        image_4 = im2uint8(imread(str2));
        image_1 = im2double(image_3);
        image_2 = im2double(image_4);

        % ����Ӱ��ռ�
        t1=clock;
        tic;
        [nonelinear_space_1,E_space_1,Max_space_1,Min_space_1,Phase_space_1]=Create_Image_space(image_1,nOctaves,Scale, ScaleValue, ratio,sigma_1,filter);
        [nonelinear_space_2,E_space_2,Max_space_2,Min_space_2,Phase_space_2]=Create_Image_space(image_2,nOctaves,Scale, ScaleValue, ratio,sigma_1,filter);
        disp(['����Ӱ��߶ȿռ仨��ʱ�䣺',num2str(toc),' ��']);

        % �������
        tic;
        [Bolb_KeyPts_1,Corner_KeyPts_1,Bolb_gradient_1,Corner_gradient_1,Bolb_angle_1,Corner_angle_1]  =  WSSF_features(nonelinear_space_1,E_space_1,Max_space_1,Min_space_1,Phase_space_1,sigma_1,ratio,Scale,nOctaves);
        [Bolb_KeyPts_2,Corner_KeyPts_2,Bolb_gradient_2,Corner_gradient_2,Bolb_angle_2,Corner_angle_2]  =  WSSF_features(nonelinear_space_2,E_space_2,Max_space_2,Min_space_2,Phase_space_2,sigma_1,ratio,Scale,nOctaves);
        disp(['��������ȡ����ʱ��:  ',num2str(toc),' ��']);

        %% ������λ��д���ļ�
        % д��kpts1���ļ�
        filenameKpts1 = fullfile(folderPath, 'WSSF_kpts1.txt');
        % �Ը���ģʽд���һ�����ݣ���֤�ļ����µĻ����ԭ������
        dlmwrite(filenameKpts1, Bolb_KeyPts_1(:,1:2), 'delimiter', ' ', 'precision', '%f');
        % ��׷��ģʽд��ڶ�������
        dlmwrite(filenameKpts1, Corner_KeyPts_1(:,1:2), 'delimiter', ' ', 'precision', '%f', '-append');
        
        % д��kpts2���ļ�
        filenameKpts2 = fullfile(folderPath, 'WSSF_kpts2.txt');
        % �Ը���ģʽд���һ������
        dlmwrite(filenameKpts2, Bolb_KeyPts_2(:,1:2), 'delimiter', ' ', 'precision', '%f');
        % ��׷��ģʽд��ڶ�������
        dlmwrite(filenameKpts2, Corner_KeyPts_2(:,1:2), 'delimiter', ' ', 'precision', '%f', '-append');


        % ��������
        tic;
        Bolb_descriptors_1 = GLOH_descriptors(Bolb_gradient_1, Bolb_angle_1, Bolb_KeyPts_1, Path_Block, ratio,sigma_1);
        Corner_descriptors_1 = GLOH_descriptors(Corner_gradient_1, Corner_angle_1, Corner_KeyPts_1, Path_Block, ratio,sigma_1);
        Bolb_descriptors_2 = GLOH_descriptors(Bolb_gradient_2, Bolb_angle_2, Bolb_KeyPts_2, Path_Block, ratio,sigma_1);
        Corner_descriptors_2 = GLOH_descriptors(Corner_gradient_2, Corner_angle_2, Corner_KeyPts_2, Path_Block, ratio,sigma_1);
        disp(['���������ӻ���ʱ��:  ',num2str(toc),' ��']); 
            
        %% ��������д���ļ�
        filenameKpts1 = fullfile(folderPath, 'WSSF_des1.txt');
        % �Ը���ģʽд���һ�����ݣ���֤�ļ����µĻ����ԭ������
        dlmwrite(filenameKpts1, Bolb_descriptors_1.des(:,1:396), 'delimiter', ' ', 'precision', '%f');
        % ��׷��ģʽд��ڶ�������
        dlmwrite(filenameKpts1, Corner_descriptors_1.des(:,1:396), 'delimiter', ' ', 'precision', '%f', '-append');
        
        % д��kpts2���ļ�
        filenameKpts2 = fullfile(folderPath, 'WSSF_des2.txt');
        % �Ը���ģʽд���һ������
        dlmwrite(filenameKpts2, Bolb_descriptors_2.des(:,1:396), 'delimiter', ' ', 'precision', '%f');
        % ��׷��ģʽд��ڶ�������
        dlmwrite(filenameKpts2, Corner_descriptors_2.des(:,1:396), 'delimiter', ' ', 'precision', '%f', '-append');


        % ����ƥ���뽻����ᴿ
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

        %% ����ƥ��ͼ��
        folderPath = 'F:\PC2-Data\UAV_GE1\Pairs_new4\results\YX\WSSF'; 
        % ����ļ����Ƿ���ڣ�����������򴴽�
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        % �����ļ�����������·��
        filename = fullfile(folderPath, sprintf('%s_WSSF.png', folderName));
        % ���浱ǰͼ�ε�ָ���ļ�
        saveas(gcf, filename);

        %% ����ƥ�������
        filename = fullfile(folderPath, 'WSSF_NCM.txt');
        fid = fopen(filename, 'a');
        fprintf(fid, '%d,%d\n', size(cleanedPoints1, 1), size(finalPoints1, 1));
        fclose(fid);
        close all;
    else
        fprintf('Files not found in %s. Skipping...\n', folderPath);
    end
end
