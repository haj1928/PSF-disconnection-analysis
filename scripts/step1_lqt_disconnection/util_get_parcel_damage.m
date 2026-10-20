# ============================================================
# 【中文备注】MATLAB 官方 bug 补丁（util_get_parcel_damage）
# 步骤: Step 1
# 说明: 本文件为研究复现代码包副本，原文件在 run62_crosssectional/lqt_patch/util_get_parcel_damage.m
#   - 修复: 官方版 L40 `lesion.img.*temp2` 整数×single 报错 → temp2=single(parcels.img)，其余逐字一致
#   - 用法: addpath 时补丁目录必须最后添加（MATLAB addpath 默认加到 path 最前，最后添加的优先级最高）
# ============================================================
%% util_get_parcel_damage —— 补丁版（修复整数类型 bug）
% =====================================================================
% 官方版 bug: L33 `temp2 = parcels.img;`（整数型）与 L40 `lesion.img.*temp2`
% （single 型）相乘报错 "整数只能与同类的整数或双精度标量值组合使用"。
% 修复: temp2 显式转 single。其余与官方逐字一致。
% 补丁位置: lqt_patch/（脚本 addpath 置于 Functions 之前 → 优先命中）
% =====================================================================

function [pcd_vect, pcd_vol] = util_get_parcel_damage(cfg)

f = filesep;

lesion = load_nii(cfg.lesion_path);
parcels = load_nii(cfg.parcel_path);

% ensure lesion is single type
lesion.img = single(lesion.img);

% Check if lesion and parcels are the same dimensions
if size(lesion.img) ~= size(parcels.img)
    error('Lesion and parcel volumes have different dimensions. Please make sure that both volumes are in the same space and have identical dimensions.');
else
    disp('Lesion and parcel volumes have identical dimensions; computing parcel damage measures.');
end

% Get parcel damage measures
pcd_vol = parcels; pcd_vol.img(:)=0; % will hold lesioned parcel values
pcd_vect = zeros(1, numel(unique(parcels.img(parcels.img>0)))); % preallocate pcd_vect
for i = 1:numel(unique(parcels.img(parcels.img > 0))) % loop through parcels
    disp(['Progress: ' num2str(i) '/' num2str(numel(unique(parcels.img(parcels.img>0)))) ' parcels evaluated']);
    temp2 = single(parcels.img); % ★ 修复: 转 single（官方版为 parcels.img 整数型）
    temp2(temp2 ~= i) = 0; % set all voxels not corresponding to this parcel to 0
    temp2(temp2 == i) = 1; % set all voxels corresponding to this parcel to 1
    pcd_vect(1,i) = (numel(find(lesion.img.*temp2))./numel(find(temp2)))*100; % get percent of parcel damaged
    pcd_vol.img(parcels.img==i) = pcd_vect(1,i); % put into volume
end
disp('Finished computing parcel damage measures');

% save parcel damage measures
cd(cfg.out_path);
if isfolder(cfg.pat_id)==0
    mkdir(cfg.pat_id);
end
cd(['.' f cfg.pat_id]);

if isfolder('Parcel_Damage') == 0
    mkdir('Parcel_Damage');
end
cd(['.' f 'Parcel_Damage']);
save_nii(pcd_vol, [cfg.pat_id '_' cfg.file_suffix '_percent_damage.nii.gz']);
save([cfg.pat_id '_' cfg.file_suffix '_percent_damage.mat'], 'pcd_vect');

end
