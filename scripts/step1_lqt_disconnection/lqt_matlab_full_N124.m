# ============================================================
# 【中文备注】断连组图计算（MATLAB 官方原版，全量 124 人）
# 步骤: Step 1
# 说明: 本文件为研究复现代码包副本，原文件在 run62_crosssectional/lqt_matlab_full_N124.m
#   - 输入: lqt_matlab_lesions_N124.csv（sub, lesion_path 两列，124 行）+ MATLAB 版 LQT 工具箱（Griffis 2021 官方原版）
#   - 输出: LQT_MATLAB_Full/{sub}/ 含 percent_parcel_SDC.mat / Disconnection_Maps/*_percent_tdi.nii.gz / Tract_Disconnection/ / Atlas/
#   - 官方来源: jcgriffis/lesion_quantification_toolkit（MATLAB 官方原版）；与 R 版同源图谱（135 区）
#   - 图谱: Schaefer2018_100Parcels_7Networks_order_plus_subcort.nii（135 区）+ HCP-842 模板
#   - 统计方法: 同 R 版（parcel 断连 = Σ断连流线/Σ图谱流线）；输出用于与 R 版逐值对比（预期 max|diff|~1e-16）
#   - 注意: 需 lqt_patch/ 补丁（官方 util_get_parcel_damage 整数×single 类型 bug）；addpath 顺序：补丁目录最后添加
#   - 运行: matlab -batch "run('F://dcmtonii//run62_crosssectional//lqt_matlab_full_N124.m')"
# ============================================================
%% MATLAB 版 LQT 全量 124 人区域级断连计算
% =====================================================================
% 2026-08-30 | MATLAB R2022b (C:/Program Files/MATLAB/R2022b)
% 工具: lesion_quantification_toolkit-main (Griffis et al. 2021 MATLAB 官方原版)
% 目的: 用 MATLAB 官方原版计算全部 124 人的区域级断连（Step 1 权威口径）
%       —— 与 5 人抽验版 (lqt_matlab_crossval.m) 同管线，扩展为全量
% 输入: lqt_matlab_lesions_N124.csv（sub, lesion_path 两列，124 行）
%       由 Python 从 FSS9_baseline_124.csv 生成（124/124 mask 已验证存在）
% 输出: LQT_MATLAB_Full/{sub}/ ... 含:
%       - percent_parcel_SDC.mat    (parcel×parcel 断连百分比矩阵)
%       - Disconnection_Maps/{sub}_Yeo7100_percent_tdi.nii.gz (断连 TDI 图)
%       - Tract_Disconnection/ ...  (70 束断连严重度)
%       - Atlas/atlas_Yeo7100_connectivity.mat
% 对照: 输出将由 Python 汇总为 124×135 矩阵，与 R 版
%       PSF_SDC_parcel_discon_N124.csv 逐值对比（预期 max|diff|~1e-16）
% 运行: 在 MATLAB 里运行本脚本 (F5)，或命令行:
%       matlab -batch "run('F:\dcmtonii\run62_crosssectional\lqt_matlab_full_N124.m')"
% =====================================================================

clc; clear;

% 固定工作目录（避免 pwd 干扰函数解析）
cd('F:\dcmtonii\run62_crosssectional');

lqt_path = 'F:\dcmtonii\run62_crosssectional\lesion_quantification_toolkit-main';
addpath([lqt_path '/Functions']);
addpath(genpath([lqt_path '/Support_Tools']));
% ★ 补丁目录必须【最后】addpath（最后添加的目录优先级最高），
%   这样 util_get_parcel_damage 命中补丁版（修整数类型 bug）而非官方版
addpath('F:\dcmtonii\run62_crosssectional\lqt_patch');

% ★ 补丁生效验证：应显示 lqt_patch\util_get_parcel_damage
fprintf('补丁验证: %s\n', which('util_get_parcel_damage'));

% ---------------- cfg 设置（与 5 人抽验版完全一致） ----------------
cfg.dsi_path = 'F:\dcmtonii\run62_crosssectional\extdata\DSI_studio\dsi_studio_64\dsi_studio.exe';
cfg.source_path = fullfile(lqt_path, 'Support_Tools', 'Tractography_Atlas', 'HCP_842');
cfg.out_path = 'F:\dcmtonii\run62_crosssectional\LQT_MATLAB_Full';
if isfolder(cfg.out_path) == 0
    mkdir(cfg.out_path);
end

% 图谱（与 R 版 100Parcels7Networks 同源，135 区）
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

% ---------------- 读 124 人列表 CSV ----------------
lesions_csv = 'F:\dcmtonii\run62_crosssectional\lqt_matlab_lesions_N124.csv';
T = readtable(lesions_csv, 'Delimiter', ',');
subs = T.sub;          % 被试 ID 列
paths = T.lesion_path; % lesion mask 路径列
n = height(T);
fprintf('读入被试数: %d\n', n);
fprintf('第 1 人: %s\n', char(subs(1)));
fprintf('最后 1 人: %s\n', char(subs(n)));

% ---------------- 全量循环 124 人 ----------------
t_start = tic;
for i = 1:n
    cfg.pat_id = char(subs(i));
    cfg.lesion_path = char(paths(i));
    fprintf('\n===== [%d/%d] %s =====\n', i, n, cfg.pat_id);

    % 1) 灰质 parcel 病灶负荷
    util_get_parcel_damage(cfg);
    % 2) 70 束断连严重度
    util_get_tract_discon(cfg);
    % 3) parcel 断连 + SSPL（内部先建 Atlas，全量共享，只建一次）
    util_get_parcel_cons(cfg);
end

fprintf('\n===== 全部 124 人完成，耗时 %.1f 分钟 =====\n', toc(t_start)/60);
fprintf('输出目录: %s\n', cfg.out_path);
