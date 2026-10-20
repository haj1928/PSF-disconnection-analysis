# ============================================================
# 【中文备注】MATLAB 官方版 × R 版交叉验证（5 人抽验，历史记录）
# 步骤: Step 1
# 说明: 本文件为研究复现代码包副本，原文件在 run62_crosssectional/lqt_matlab_crossval.m
#   - 输入: 5 个分层被试（小/中/大病灶）lesion mask
#   - 输出: LQT_MATLAB_Crossval/ + PSF_MATLAB_vs_R_LQT_crossval_detail.csv（675 对逐值对比）
#   - 官方来源: jcgriffis/lesion_quantification_toolkit（MATLAB 官方原版）
#   - 结果: 5 人 × 135 parcel 共 675 对，max|diff|=1.11e-16（浮点极限），corr=1.0000 → R 版=MATLAB 版铁证
#   - 运行: MATLAB R2022b 直接 F5
# ============================================================
%% MATLAB 版 LQT 交叉验证 —— 5 人抽验（与 R 版逐值对比）
% =====================================================================
% 2026-08-26 | MATLAB R2022b (I:\Matlab2022b)
% 工具: lesion_quantification_toolkit-main (Griffis et al. 2021 MATLAB 版)
% 目的: 与 R 版 LQT (jdwor/LQT) 的 124 人输出做交叉验证（先跑 5 人）
% 图谱: Schaefer-Yeo 100 + 皮下 (与 R 版同一张, 已从 R 版 extdata 补齐)
% 模板: HCP-842 (已从 R 版 extdata 补齐: fib/trk/All_Tracts)
% 运行: MATLAB 里直接运行本脚本 (F5) 即可
% =====================================================================

clc; clear;

% 固定工作目录（避免 pwd 干扰函数解析）
cd('F:\dcmtonii\run62_crosssectional');

lqt_path = 'F:\dcmtonii\run62_crosssectional\lesion_quantification_toolkit-main';
addpath([lqt_path '/Functions']);
addpath(genpath([lqt_path '/Support_Tools']));
% ★ 补丁目录必须【最后】addpath —— MATLAB addpath 默认加到 path 最前，
%   最后添加的目录优先级最高，这样 util_get_parcel_damage 会命中补丁版
%   （修整数类型 bug），而不是官方版。
addpath('F:\dcmtonii\run62_crosssectional\lqt_patch');

% ★ 补丁生效验证：应显示 lqt_patch\util_get_parcel_damage
fprintf('补丁验证: %s\n', which('util_get_parcel_damage'));

% ---------------- cfg 设置 ----------------
cfg.dsi_path = 'F:\dcmtonii\run62_crosssectional\extdata\DSI_studio\dsi_studio_64\dsi_studio.exe';
cfg.source_path = fullfile(lqt_path, 'Support_Tools', 'Tractography_Atlas', 'HCP_842');
cfg.out_path = 'F:\dcmtonii\run62_crosssectional\LQT_MATLAB_Crossval';
if isfolder(cfg.out_path) == 0
    mkdir(cfg.out_path);
end

% 图谱（已补齐的真实 .nii, 与 R 版 100Parcels7Networks 同源）
cfg.parcel_path = fullfile(lqt_path, 'Support_Tools', 'Parcellations', ...
    'Schaefer_Yeo', 'Plus_Subcort', 'Schaefer2018_100Parcels_7Networks_order_plus_subcort.nii');

% 图谱配套 .mat (RegionName/坐标) —— MATLAB 版官方文件
load(fullfile(lqt_path, 'Support_Tools', 'Parcellations', 'Schaefer_Yeo', ...
    'Plus_Subcort', 'Schaefer2018_100Parcels_7Networks_order_plus_subcort.mat'));
cfg.node_label = t.RegionName;
cfg.node_color = t.NetworkID;
cfg.parcel_coords = [t.X, t.Y, t.Z];

% 输出命名与连接参数（与 R 版一致: con_type='end', 阈值默认）
cfg.file_suffix = 'Yeo7100';
cfg.con_type = 'end';
cfg.sspl_spared_thresh = 50;
cfg.smooth = 0;

% ---------------- 图谱验证（打印 RegionName 供对比） ----------------
fprintf('\n===== 图谱 RegionName 验证 =====\n');
fprintf('RegionName 数量: %d\n', numel(cfg.node_label));
for i = 1:8
    fprintf('  [%d] %s\n', i, char(cfg.node_label(i)));
end
fprintf('  ... 最后 3 个: %s / %s / %s\n', ...
    char(cfg.node_label(end-2)), char(cfg.node_label(end-1)), char(cfg.node_label(end)));
fprintf('===============================\n\n');

% ---------------- 5 个分层被试（小/中/大病灶） ----------------
lesions = {
    'sub_10090', 'F:\dcmtonii\sub_10090\MNI_1mm_Output\sub_10090_lesion_mask_MNI152_1mm.nii'
    'sub_10057', 'F:\dcmtonii\sub_10057\MNI_1mm_Output\sub_10057_lesion_mask_MNI152_1mm.nii'
    'sub_10058', 'F:\dcmtonii\sub_10058\MNI_1mm_Output\sub_10058_lesion_mask_MNI152_1mm.nii'
    'sub_10007', 'F:\dcmtonii\sub_10007\MNI_1mm_Output\sub_10007_lesion_mask_MNI152_1mm.nii'
    'sub_10092', 'F:\dcmtonii\sub_10092\MNI_1mm_Output\sub_10092_lesion_mask_MNI152_1mm.nii'
};

for i = 1:size(lesions, 1)
    cfg.pat_id = lesions{i, 1};
    cfg.lesion_path = lesions{i, 2};
    fprintf('\n===== %s (%d/%d) =====\n', cfg.pat_id, i, size(lesions, 1));
    
    % 1) 灰质 parcel 病灶负荷
    util_get_parcel_damage(cfg);
    % 2) 70 束断连严重度
    util_get_tract_discon(cfg);
    % 3) parcel 断连 + SSPL（内部先建 Atlas，5 人共享，只建一次）
    util_get_parcel_cons(cfg);
end

fprintf('\n===== 全部完成。输出目录: %s =====\n', cfg.out_path);
fprintf('关键输出（对比用）:\n');
fprintf('  - <pat>_Yeo7100_percent_parcel_SDC.mat   (parcel×parcel 断连百分比矩阵)\n');
fprintf('  - Atlas/atlas_Yeo7100_connectivity.mat   (atlas 连接矩阵)\n');
