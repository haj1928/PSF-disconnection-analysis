# =============================================================================
# psf_sdc_lqt_n124_v2.R —— 断连组图 124 人重跑（修复版 v2.4）
# =============================================================================
# 【修复 1】跳过 get_parcel_discon 的 "Save TDI map" 段
#   ★ 根因（已读 LQT 源码 + oro.nifti 文档确认，非猜测）：
#     LQT::get_parcel_discon 在 "Save TDI map" 阶段用 neurobase::readnii() 读取
#     DSI Studio 生成的 .tdi.nii.gz。该文件 NIfTI header 的 magic 字段为 "ni1"
#     （dual-file 布局标记），而 readnii() 底层 oro.nifti::readNIfTI 默认按
#     single-file 读取，命中检查并抛 "This is not a one-file NIfTI format"。
#     （oro.nifti 文档实证：readNIfTI 读 magic="ni1" 的文件即报此错）
#   ★ 为什么可以安全跳过：
#     TDI 仅用于生成 percent_tdi.nii.gz 可视化图。我们需要的 5 个 CSV
#     （parcel_discon / parcel_damage / tract_discon / net_discon / parc2parc）
#     均不依赖 TDI——已核对 LQT 源码：
#       - parcel_overall_discon.csv 在 TDI 段之前写出（get_parcel_discon.R L150）
#       - get_patient_sspl.R 只读 _percent_parcel_mats.RData（不读 TDI）
#       - compile_data.R 只汇总 CSV（不读 TDI）
# 【修复 2】parcel_path 必须用 LQT 内置 extdata 路径（.nii.gz）
#   ★ 根因（已读 LQT 源码确认）：
#     create_cfg_object.R L80-94：当 parcel_path 是自定义路径（不含 "extdata"）
#     时进入第一个分支，但该分支【漏写了 t=read.csv(...)】，导致 t 仍是 R 内置
#     转置函数（closure），随后 t$RegionName 报 "类别为'closure'的对象不可以取子集"。
#     而 extdata 分支（L95-100）和默认分支（L106-114）都有 read.csv，正常。
#     → 因此 parcel_path 一律用 system.file 内置路径，不要用自定义路径。
# 【修复 3】compile_data 用全局版（mclapply → 串行 lapply）
#   ★ 根因（已读 LQT 源码确认）：
#     compile_data.R L33/38/43/45/50/55/60 共 6 处用 parallel::mclapply，
#     Windows 不支持 fork 并行 → "Windows不支持'mc.cores' > 1"。
#     汇总阶段数据量小，改串行 lapply 完全等价且最稳。
# 【修复 4】汇总前自动检查输出完整性并补跑缺失被试
#   ★ 背景：pbmclapply 对失败被试不中断（进度 100% ≠ 全部成功），个别被试可能
#     缺 Parcel_Damage csv / Tract_Disconnection csv / Parcel_Disconnection RData，
#     导致 compile_data 内 read.csv(character(0)) 报 "description 参数无效"。
#     脚本在第 5 节前自动扫描 124 人三件套，缺失则只对缺失被试补跑（最多 3 轮）。
# 【运行方式】与 v1 相同：RStudio Source 或 Rscript。
# 【时间节省】本次已跑完 parcel damage (~52min) + tract discon (~1h54min) +
#   atlas 创建 + atlas SSPL + get_parcel_discon (~10min) + get_patient_sspl
#   （约 4 小时），本脚本自动检测并跳过这些已完成的步骤，从崩点继续。
# =============================================================================

suppressMessages({
  library(LQT)
  # ★ 补丁版 get_parcel_discon 定义在全局环境，其内部调用的函数需显式加载：
  #   pbmclapply → pbmcapply；readMat → R.matlab；readnii → neurobase
  #   （LQT 原版函数在包 namespace 内无需此步，但全局补丁函数必须显式 library）
  library(pbmcapply)
  library(R.matlab)
  library(neurobase)
  library(dplyr)
  library(readr)
})

# ---------------- 0. 路径配置（★ 运行前确认） ----------------
mask_root   <- "F:/psf_rerun/PSF_lesion_masks_N124"   # 124 个 mask 所在根目录
out_root    <- "F:/psf_rerun"                          # 输出目录
# parcel_path 必须用 LQT 内置 extdata 路径（.nii.gz）——不要改成自定义路径！
# 原因：create_cfg_object.R 对"自定义路径"分支存在 bug（漏 read.csv，t 未赋值
# → t$RegionName 报 closure 错）；extdata 分支正常。且 DSI Studio 读 .nii.gz
# 已验证可行（07-11 65 人版与本次 atlas 创建均用 .nii.gz 成功）。
parcel_path <- system.file("extdata", "Schaefer_Yeo_Plus_Subcort",
                           "100Parcels7Networks.nii.gz", package = "LQT")

out_path <- file.path(out_root, "LQT_SDC_Results_N124")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# DSI Studio：留 NULL 让 LQT 使用内置（extdata/DSI_studio）
dsi_path <- NULL

# ---------------- 1. 列出 124 个病灶 mask ----------------
lesion_paths <- list.files(mask_root,
                           pattern = "lesion_mask_MNI152_1mm\\.nii$",
                           recursive = TRUE, full.names = TRUE)
pat_ids <- basename(dirname(lesion_paths))     # 每个 mask 的父目录名 = sub_xxx

if (length(pat_ids) != 124) {
  stop("mask 数量=", length(pat_ids), "，应为 124。请检查 mask_root 是否完整解压。")
}
dup <- pat_ids[duplicated(pat_ids)]
if (length(dup) > 0) stop("存在重复被试: ", paste(dup, collapse = ", "))
message("找到 ", length(pat_ids), " 个病灶 mask（124 人全覆盖）")

# ---------------- 2. 构建 cfg（与原版完全一致） ----------------
cfg <- create_cfg_object(pat_ids     = pat_ids,
                         lesion_paths = lesion_paths,
                         out_path     = out_path,
                         parcel_path  = parcel_path,
                         dsi_path     = dsi_path,
                         con_type     = "end")   # 'end' = 基于纤维束端点（保守）

# ---------------- 3. 自动诊断（验证 TDI 根因，不中断） ----------------
acfg <- cfg[[1]]
diag_pd <- file.path(out_path, acfg$pat_id, "Parcel_Disconnection")
cat("\n=== 诊断：DSI 生成的 TDI 文件能否被 readnii 读取 ===\n")
if (dir.exists(diag_pd)) {
  tdi_files <- list.files(diag_pd, pattern = "\\.tdi\\.nii\\.gz$", full.names = TRUE)
  cat("  ", acfg$pat_id, "Parcel_Disconnection 中 .tdi.nii.gz:", length(tdi_files), "个\n")
  if (length(tdi_files) > 0) {
    res <- tryCatch({
      img <- readnii(tdi_files[1])
      paste0("readnii OK, magic=", img@magic)
    }, error = function(e) paste0("readnii 报错: ", e$message))
    cat("  -> ", res, "\n")
    if (grepl("one-file NIfTI", res)) {
      cat("  -> 确认根因：TDI header magic=ni1，跳过 TDI 块即可修复 ✓\n")
    }
  } else {
    cat("  -> 该目录无 .tdi.nii.gz（上次未跑到此步），跳过\n")
  }
} else {
  cat("  -> 该被试目录尚不存在（首次进入 get_parcel_discon），跳过\n")
}
cat("====================================================\n\n")

# ---------------- ★ 补丁：覆盖 LQT 的 get_parcel_discon（跳过 TDI 块） ----------------
# 与原版 LQT::get_parcel_discon 逐字一致，仅将 "Save TDI map" 段替换为"清理+跳过"。
# 定义为全局函数，遮蔽包内同名函数；由于我们不调用 get_parcel_cons（其内部
# 会解析到包内版本），而是手动拆分调用，因此实际执行的是本补丁版。
get_parcel_discon <- function(cfg, cores = 1, verbose = T) {
  if (is.null(cfg$pat_id)) {

    at.path = paste0(cfg[[1]]$out_path, "/Atlas")
    if (!dir.exists(at.path)) {
      stop('Atlas folder does not exist in output directory; cannot convert streamline counts or TDI values to % disconnection. Please run get_parcel_atlas before attempting to create patient disconnection measures.');
    }
    cat('Computing parcel disconnection measures\n')
    out = pbmclapply(cfg, get_parcel_discon, verbose = F, mc.cores = cores)
    cat('Finished computing parcel disconnection measures.\n')

  } else {

    at.path = paste0(cfg$out_path, "/Atlas")
    if (!dir.exists(at.path)) {
      stop('Atlas folder does not exist in output directory; cannot convert streamline counts or TDI values to % disconnection. Please run get_parcel_atlas before attempting to create patient disconnection measures.');
    }
    pat.path = paste0(cfg$out_path, "/", cfg$pat_id)
    if (!dir.exists(pat.path)) {
      dir.create(pat.path)
    }
    pd.path = paste0(pat.path, "/Parcel_Disconnection")
    dm.path = paste0(pat.path, "/Disconnection_Maps")
    if (!dir.exists(pd.path)) {
      dir.create(pd.path)
    }
    if (!dir.exists(dm.path)) {
      dir.create(dm.path)
    }

    out_file = pd.path

    # create disconnection matrix, .trk file, and TDI map
    out = tryCatch({
      suppressWarnings(system(paste0('! ', cfg$dsi_path, ' --action=ana --source=', cfg$source_path, '/HCP842_1mm.fib.gz',
                                     ' --tract=', cfg$source_path, '/all_tracts_1mm.trk.gz', ' --roi=', cfg$lesion_path,
                                     ' --output=', out_file, ' --connectivity=', cfg$parcel_path, ' --connectivity_type=',
                                     cfg$con_type, ' --connectivity_threshold=0 --export=tdi'), intern = T))
    }, error = function(e) {
      if (grepl(":", cfg$lesion_path)) {
        drive = paste0(strsplit(cfg$lesion_path, ":")[[1]][1], ":")
        fsplit = strsplit(cfg$lesion_path, "/")[[1]]
        direc = paste(fsplit[-length(fsplit)], collapse = "/")
        lp = tail(fsplit, 1)
        suppressWarnings(shell(paste0(drive, " &  cd ", direc, " & ", cfg$dsi_path, ' --action=ana --source=', cfg$source_path,
                                      '/HCP842_1mm.fib.gz', ' --tract=', cfg$source_path, '/all_tracts_1mm.trk.gz', ' --roi=',
                                      lp, ' --output=', out_file, ' --connectivity=', cfg$parcel_path, ' --connectivity_type=',
                                      cfg$con_type, ' --connectivity_threshold=0 --export=tdi'), intern = T))
      } else {
        suppressWarnings(shell(paste0(cfg$dsi_path, ' --action=ana --source=', cfg$source_path, '/HCP842_1mm.fib.gz',
                                      ' --tract=', cfg$source_path, '/all_tracts_1mm.trk.gz', ' --roi=', cfg$lesion_path,
                                      ' --output=', out_file, ' --connectivity=', cfg$parcel_path, ' --connectivity_type=',
                                      cfg$con_type, ' --connectivity_threshold=0 --export=tdi'), intern = T))
      }
    })

    oldnames = list.files(pd.path, pattern = "all_tracts_1mm", full.names = T)
    newnames = gsub("all_tracts_1mm", paste0(cfg$pat_id, "_", cfg$file_suffix), oldnames)
    file.rename(oldnames, newnames)

    # resave connectivity file (raw streamline counts)
    matfile = paste0(pd.path, "/", list.files(pd.path, pattern = "connectivity\\.mat$"))
    mat = readMat(matfile)
    disconnectivity = mat$connectivity
    rownames(disconnectivity) = cfg$node_label
    colnames(disconnectivity) = cfg$node_label
    node_label = cfg$node_label
    node_group = cfg$node_group
    file.remove(matfile); rm(mat)

    # re-save network measures file (raw streamline counts)
    netfile = paste0(pd.path, "/", list.files(pd.path, pattern = "network_measures\\.txt$"))
    global = read.table(netfile, sep = "\t")[1:27, ]
    colnames(global) = c("Measure", "Value")
    local = read.table(netfile, sep = "\t", skip = 27, header = T)[, -(length(node_label) + 2)]
    colnames(local)[1] = "Measure"
    colnames(local)[-1] = cfg$node_label
    file.remove(netfile)

    save(global, local, file = paste0(pd.path, "/", cfg$pat_id, "_",
                                      cfg$file_suffix, "_network_measures.RData"))
    save(disconnectivity, node_label, node_group,
         file = paste0(pd.path, "/", cfg$pat_id, "_",
                       cfg$file_suffix, "_disconnectivity.RData"))

    # load atlas SC matrix
    pat_con = disconnectivity; rm(disconnectivity)
    load(paste0(at.path, '/atlas_', cfg$file_suffix, '_connectivity.RData'))
    atlas_con = connectivity; rm(connectivity)

    # convert patient matrix to % disconnection and save
    pct_sdc_matrix = 100 * (pat_con / atlas_con)
    pct_sdc_matrix[is.na(pct_sdc_matrix)] = 0
    pct_spared_sc_matrix = 100 * ((atlas_con - pat_con) / atlas_con)
    pct_spared_sc_matrix[is.na(pct_spared_sc_matrix)] = 0
    pct_spared_sc_matrix[is.infinite(pct_spared_sc_matrix)] = 0

    rownames(pct_sdc_matrix) = cfg$node_label
    rownames(pct_sdc_matrix) = cfg$node_label
    colnames(pct_spared_sc_matrix) = cfg$node_label
    colnames(pct_spared_sc_matrix) = cfg$node_label

    fname = paste0(pd.path, "/", cfg$pat_id, "_", cfg$file_suffix, "_percent_parcel_mats.RData")
    save(pct_sdc_matrix, pct_spared_sc_matrix, file = fname)

    ### save parcel-level disconnection percentages

    parc.discon = data.frame(Parcel = cfg$node_label, ParcelGroup = cfg$node_group)
    parc.discon$Parcel = gsub(" \\(Right\\)", "_R", parc.discon$Parcel)
    parc.discon$Parcel = gsub(" \\(Left\\)", "_L", parc.discon$Parcel)
    parc.discon$Parcel = gsub("Lenticular nucleus, p",
                              "Lent_Nuc_P", parc.discon$Parcel)
    parc.discon$Parcel = gsub("RH_", "R_", parc.discon$Parcel)
    parc.discon$Parcel = gsub("LH_", "L_", parc.discon$Parcel)

    pd.path = paste0(cfg$out_path, "/", cfg$pat_id, "/Parcel_Disconnection")
    load(paste0(pd.path, "/", cfg$pat_id, "_", cfg$file_suffix,
                "_percent_parcel_mats.RData"))
    at.path = paste0(cfg$out_path, "/Atlas")
    load(paste0(at.path, '/atlas_', cfg$file_suffix, '_connectivity.RData'))

    at_con = connectivity; rm(connectivity)
    pat_sdc = at_con * pct_sdc_matrix / 100

    parc.discon$ParcelDiscon = apply(pat_sdc, 1, sum) / apply(at_con, 1, sum)
    write.csv(parc.discon, paste0(pd.path, "/", cfg$pat_id, "_",
                                  cfg$file_suffix, "_parcel_overall_discon.csv"))

    ### output .node and .edge files for external viewers (e.g. MRIcroGL)

    # get parcel coordinates if not supplied
    if (is.null(cfg$parcel_coords)) {
      cat('Warning: Parcel coordinates not provided, attempting to estimate coordinates from parcellation image\n');
      cfg$parcel_coords = util_get_coords(cfg)
    }
    NodePos = cfg$parcel_coords

    # write out .edge files
    write(round(t(pct_sdc_matrix), 4),
          paste0(pd.path, "/", cfg$pat_id, "_", cfg$file_suffix,
                 "_percent_parcel_SDC.edge"),
          ncolumns = ncol(pct_sdc_matrix), sep = "\t")
    write(round(t(pct_spared_sc_matrix), 4),
          paste0(pd.path, "/", cfg$pat_id, "_", cfg$file_suffix,
                 "_percent_parcel_spared_SC.edge"),
          ncolumns = ncol(pct_spared_sc_matrix), sep = "\t")

    # write out .node files
    NodeSize = apply(pct_sdc_matrix, 2, sum, na.rm = T) / apply((atlas_con > 0) * 100, 2, sum, na.rm = T) # size nodes according to % maximum # of affected connections
    NodeSize[is.na(NodeSize)] = 0
    NodeSize2 = apply(pct_spared_sc_matrix, 2, sum, na.rm = T) / apply((atlas_con > 0) * 100, 2, sum, na.rm = T)
    NodeSize2[is.na(NodeSize2)] = 0
    if (is.null(cfg$node_color)) {
      NodeColor = NodeSize
    } else if (length(cfg$node_color) != length(NodeSize)) {
      NodeColor = NodeSize
    } else {
      NodeColor = cfg$node_color
    }
    if (is.null(cfg$node_label)) {
      NodeLabel = as.character(1:nrow(pct_sdc_matrix))
    } else if (length(cfg$node_label) != length(NodeSize)) {
      NodeLabel = as.character(1:nrow(pct_sdc_matrix))
    } else {
      NodeLabel = cfg$node_label
    }

    # construct output filename and save
    OutputFile = paste0(pd.path, "/", cfg$pat_id, "_", cfg$file_suffix, "_percent_parcel_SDC.node")
    C = cbind(NodePos, NodeColor, NodeSize, NodeLabel)
    write(t(C), OutputFile, sep = "\t")

    OutputFile = paste0(pd.path, "/", cfg$pat_id, "_", cfg$file_suffix, "_percent_parcel_spared_SC.node")
    C = cbind(NodePos, NodeColor, NodeSize2, NodeLabel)
    write(t(C), OutputFile, sep = "\t")

    ### Save TDI map —— ★ 补丁：跳过（原版在此崩）
    # 根因：DSI Studio 生成的 .tdi.nii.gz 的 NIfTI header magic="ni1"（dual-file 标记），
    # oro.nifti::readNIfTI 默认 single-file 读取 → "This is not a one-file NIfTI format"。
    # TDI 仅用于生成 percent_tdi.nii.gz 可视化图，5 个 CSV 输出均不依赖（已核对
    # compile_data.R / get_patient_sspl.R 源码），故跳过并清理 DSI 中间文件。
    trk_file = list.files(pd.path, pattern = "\\.trk\\.gz$")
    if (length(trk_file) > 0) { file.remove(paste0(pd.path, "/", trk_file)) }
    tdi_files = list.files(pd.path, pattern = "\\.tdi\\.nii\\.gz$")
    if (length(tdi_files) > 0) { file.remove(paste0(pd.path, "/", tdi_files)) }

    if (verbose == T) { cat('Finished computing parcel disconnection measures.\n') }
  }
}

# ---------------- ★ 补丁：覆盖 LQT 的 compile_data（mclapply → 串行 lapply） ----------------
# 原版 compile_data.R 用 parallel::mclapply（L33/38/43/45/50/55/60），Windows 不支持
# fork 并行 → "Windows不支持'mc.cores' > 1"。汇总阶段数据量小（124 人），串行 lapply
# 完全等价且最稳。内部调用的 pdam/tdis/p2pnet/ndis/pdis/n2ndis/p2pdis 是 LQT 未导出
# 函数，用 LQT::: 访问。其余逻辑与源码逐字一致。
compile_data <- function(cfg, cores = 1) {
  pat_ids = unlist(lapply(cfg, `[[`, 1))
  out_path = cfg[[1]]$out_path
  suffix = cfg[[1]]$file_suffix
  at.path = file.path(out_path, "Atlas")
  pat.paths = file.path(out_path, pat_ids)
  ID = pat_ids
  load(paste0(at.path, '/atlas_', suffix, '_connectivity.RData'))
  atlas = connectivity; rm(connectivity)

  pdam_f    <- LQT:::pdam
  tdis_f    <- LQT:::tdis
  p2pnet_f  <- LQT:::p2pnet
  ndis_f    <- LQT:::ndis
  pdis_f    <- LQT:::pdis
  n2ndis_f  <- LQT:::n2ndis
  p2pdis_f  <- LQT:::p2pdis

  parc.damage = lapply(pat.paths, pdam_f, suffix)
  parc.damage = do.call(rbind, parc.damage)
  parc.damage = as.data.frame(parc.damage)
  parc.damage = cbind(ID, parc.damage)

  tract.discon = lapply(pat.paths, tdis_f)
  tract.discon = do.call(rbind, tract.discon)
  tract.discon = as.data.frame(tract.discon)
  tract.discon = cbind(ID, tract.discon)

  patnets = lapply(pat.paths, p2pnet_f, suffix)

  net.discon = lapply(patnets, ndis_f, atlas, node_group)
  net.discon = do.call(rbind, net.discon)
  net.discon = as.data.frame(net.discon)
  net.discon = cbind(ID, net.discon)

  parc.discon = lapply(patnets, pdis_f, atlas)
  parc.discon = do.call(rbind, parc.discon)
  parc.discon = as.data.frame(parc.discon)
  parc.discon = cbind(ID, parc.discon)

  net2net.discon = lapply(patnets, n2ndis_f, atlas, node_group)
  net2net.discon = do.call(rbind, net2net.discon)
  net2net.discon = as.data.frame(net2net.discon)
  net2net.discon = cbind(ID, net2net.discon)

  parc2parc.discon = lapply(patnets, p2pdis_f, atlas)
  parc2parc.discon = do.call(rbind, parc2parc.discon)
  allnas = apply(parc2parc.discon, 2, function(x) mean(is.na(x)) == 1)
  parc2parc.discon = as.data.frame(parc2parc.discon[, !allnas])
  parc2parc.discon = cbind(ID, parc2parc.discon)

  output = list(parc.damage = parc.damage, tract.discon = tract.discon,
                net.discon = net.discon, parc.discon = parc.discon,
                net2net.discon = net2net.discon,
                parc2parc.discon = parc2parc.discon)
  return(output)
}

# ---------------- 4. 计算（跳过已完成步骤，从崩点继续） ----------------
# ✅ 本次已跑完（跳过，避免重复计算约 4 小时）：
#   get_parcel_damage(cfg, cores=2)   # 已完成 ~52min
#   get_tract_discon(cfg, cores=2)    # 已完成 ~1h54min
#   get_parcel_atlas / get_atlas_sspl # 已完成（Atlas 目录存在）
# 以下检查确保 atlas 若缺失则自动补建：
at.path <- file.path(out_path, "Atlas")
if (!dir.exists(at.path)) {
  message("Atlas 目录不存在，创建 atlas 连接矩阵与 SSPL...")
  get_parcel_atlas(cfg)
  get_atlas_sspl(cfg)
} else {
  con.file <- list.files(at.path, pattern = paste0(acfg$file_suffix, "_connectivity\\.RData$"))
  if (length(con.file) == 0) {
    message("Atlas connectivity 缺失，重建...")
    get_parcel_atlas(cfg)
  }
  sspl.file <- list.files(at.path, pattern = paste0(acfg$file_suffix, "_SSPL_matrix\\.RData$"))
  if (length(sspl.file) == 0) {
    message("Atlas SSPL 缺失，重建...")
    get_atlas_sspl(cfg)
  }
}

# ---- get_parcel_discon（补丁版，跳过 TDI）：若已完成则跳过重跑 ----
first_dir <- file.path(out_path, pat_ids[1])
first_pd  <- file.path(first_dir, "Parcel_Disconnection")
discon_done <- dir.exists(first_pd) &&
  length(list.files(first_pd, pattern = "_percent_parcel_mats\\.RData$")) > 0
if (discon_done) {
  message("检测到 get_parcel_discon 已完成（", length(pat_ids), " 人），跳过重跑")
} else {
  # 仅当未完成时，清理首被试可能的半成品再干净重跑
  if (dir.exists(first_dir)) {
    message("清理上次失败残留: ", first_dir)
    unlink(first_dir, recursive = TRUE)
  }
  message("开始 get_parcel_discon（124 人，预计 1-2 小时）...")
  get_parcel_discon(cfg, cores = 2)
}

# ---- get_patient_sspl（LQT 原版）：若已完成则跳过 ----
first_ps <- file.path(first_dir, "Parcel_SSPL")
sspl_done <- dir.exists(first_ps) &&
  length(list.files(first_ps, pattern = "_SSPL_matrices\\.RData$")) > 0
if (sspl_done) {
  message("检测到 get_patient_sspl 已完成，跳过")
} else {
  message("开始 get_patient_sspl（124 人）...")
  get_patient_sspl(cfg, cores = 2)
}

# ---------------- 4.5 输出完整性检查 + 自动补跑缺失被试 ----------------
# 背景：pbmclapply 对失败被试不中断（进度 100% ≠ 全部成功），个别被试可能缺文件，
# 导致 compile_data 里 read.csv(character(0)) 报 "description 参数无效"。
# 这里在汇总前扫描每被试三件套，缺失则自动只补跑该被试，最多 3 轮。
suffix_cfg <- cfg[[1]]$file_suffix
check_completeness <- function() {
  out <- data.frame(sub = pat_ids, damage = NA_integer_, tract = NA_integer_, discon = NA_integer_)
  for (i in seq_along(pat_ids)) {
    s <- pat_ids[i]
    pd <- file.path(out_path, s, "Parcel_Damage")
    td <- file.path(out_path, s, "Tract_Disconnection")
    pc <- file.path(out_path, s, "Parcel_Disconnection")
    out$damage[i] <- if (dir.exists(pd)) length(list.files(pd, pattern = paste0(suffix_cfg, "_percent_damage\\.csv$"))) else 0L
    out$tract[i]  <- if (dir.exists(td)) length(list.files(td, pattern = "\\.csv$")) else 0L
    out$discon[i] <- if (dir.exists(pc)) length(list.files(pc, pattern = paste0(suffix_cfg, "_disconnectivity\\.RData$"))) else 0L
  }
  out
}

for (rnd in 1:3) {
  comp <- check_completeness()
  miss <- which(comp$damage == 0 | comp$tract == 0 | comp$discon == 0)
  if (length(miss) == 0) break
  message("【完整性检查】第 ", rnd, " 轮：", length(miss), " 个被试缺文件，自动补跑中...")
  print(comp[miss, ])
  miss_damage <- which(comp$damage == 0)
  miss_tract  <- which(comp$tract == 0)
  miss_discon <- which(comp$discon == 0)
  if (length(miss_damage) > 0) { message("补跑 get_parcel_damage: ", length(miss_damage), " 人"); get_parcel_damage(cfg[miss_damage], cores = 2) }
  if (length(miss_tract)  > 0) { message("补跑 get_tract_discon: ", length(miss_tract),  " 人"); get_tract_discon(cfg[miss_tract], cores = 2) }
  if (length(miss_discon) > 0) { message("补跑 get_parcel_discon: ", length(miss_discon), " 人"); get_parcel_discon(cfg[miss_discon], cores = 2) }
}
comp <- check_completeness()
miss <- which(comp$damage == 0 | comp$tract == 0 | comp$discon == 0)
if (length(miss) > 0) {
  message("!! 仍有 ", length(miss), " 个被试文件缺失（补跑也失败），无法汇总：")
  print(comp[miss, ])
  stop("输出不完整，请检查缺失被试的 DSI Studio 报错。")
}
message("输出完整性检查通过：", nrow(comp), " 人三件套齐全")

# ---------------- 5. 汇总并写出 _N124.csv ----------------
message("汇总并写出 _N124.csv ...")
data <- compile_data(cfg, cores = 2)

write_csv(data$parc.discon,      file.path(out_root, "PSF_SDC_parcel_discon_N124.csv"))
write_csv(data$parc.damage,      file.path(out_root, "PSF_SDC_parcel_damage_N124.csv"))
write_csv(data$tract.discon,     file.path(out_root, "PSF_SDC_tract_discon_N124.csv"))
write_csv(data$net.discon,       file.path(out_root, "PSF_SDC_net_discon_N124.csv"))
write_csv(data$parc2parc.discon, file.path(out_root, "PSF_SDC_parc2parc_discon_N124.csv"))

message("完成。SDC 特征（_N124）已写入 ", out_root)
message("主特征矩阵：PSF_SDC_parcel_discon_N124.csv（124 人 × parcel）")
