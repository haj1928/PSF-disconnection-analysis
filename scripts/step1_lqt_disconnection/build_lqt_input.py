# -*- coding: utf-8 -*-
"""
code_build_lqt_input_newmasks_20260915.py

为"新掩膜"下游重算建立 LQT 输入清单，并区分"必须重跑"与"可复用既有 LQT 输出"。

背景
----
- 旧 LQT 输入 lqt_matlab_lesions_N124.csv 指向 sub_*/MNI_1mm_Output/*_MNI152_1mm.nii，
  其 mtime 呈两批：2026-07-09（63 例，厚层 T1 路线）与 2026-08-07（61 例，薄层路线）。
- 2026-09-14 的重跑只改了那 63 例（新掩膜在 sub_*/MNI/lesion_mni_*.nii，interp=0 薄层路线；
  sub_10005 为其自身经验证正确的厚层路线结果）。
- 61 例（8 月薄层批次）的掩膜没有变化 → 可复用既有 LQT_MATLAB_Full/<sub>/ 输出，无需重跑。

输出
----
- lqt_matlab_lesions_N124_newmasks.csv   : sub,lesion_path（统一指向"当前口径"掩膜）
- mask_switch_plan_20260915.csv          : 逐例记录（新旧路径、批次、Dice、是否复用 LQT）
"""
import io
import os
import numpy as np
import nibabel as nib
import pandas as pd
from nibabel.processing import resample_from_to

R = r'F:/dcmtonii'
RUN = os.path.join(R, 'run62_crosssectional')
OUT_DIR = os.path.join(RUN, 'PSF_pls_behavior', '断连稿写作_20260907', '修订_revision_20260910')
CSV_IN = os.path.join(RUN, 'lqt_matlab_lesions_N124.csv')
CSV_OUT = os.path.join(RUN, 'lqt_matlab_lesions_N124_newmasks.csv')
PLAN = os.path.join(OUT_DIR, 'mask_switch_plan_20260915.csv')

d = pd.read_csv(CSV_IN)
assert list(d.columns) == ['sub', 'lesion_path'] and len(d) == 124


def load_bin(path):
    im = nib.load(path)
    a = np.asanyarray(im.dataobj)
    return im, (a > 0.5)


def dice_grid(a_bin, a_img, b_bin, b_img):
    """把 b 重采样到 a 的网格后算 Dice（网格不同也没问题）"""
    if a_bin.shape == b_bin.shape and np.allclose(a_img.affine, b_img.affine):
        bb = b_bin
    else:
        bb = np.asanyarray(resample_from_to(
            nib.Nifti1Image(b_bin.astype(np.uint8), b_img.affine),
            (a_img.shape, a_img.affine), order=0).dataobj) > 0.5
    if a_bin.sum() + bb.sum() == 0:
        return float('nan')
    return 2.0 * (a_bin & bb).sum() / (a_bin.sum() + bb.sum())


rows = []
for r in d.itertuples():
    s = str(r.sub)
    old = str(r.lesion_path)
    new = os.path.join(R, s, 'MNI', 'lesion_mni_%s.nii' % s)
    use_new = os.path.isfile(new)
    use = new if use_new else old
    old_dt = pd.Timestamp(os.path.getmtime(old), unit='s')
    batch = '2026-07-09_thick63' if old_dt.strftime('%Y-%m-%d') == '2026-07-09' else '2026-08-07_thin61'
    reuse = batch == '2026-08-07_thin61'          # 61 例掩膜未变 → 复用既有 LQT 输出
    dm = float('nan')
    if use_new:
        try:
            oi, ob = load_bin(old)
            ni, nb = load_bin(new)
            dm = dice_grid(nb, ni, ob, oi)
        except Exception:
            dm = float('nan')
    rows.append(dict(subject=s, batch=batch, reuse_lqt=reuse,
                     mask_used=use, new_path_exists=use_new, old_path=old,
                     dice_new_vs_old=dm,
                     old_mtime=old_dt.strftime('%Y-%m-%d %H:%M')))

plan = pd.DataFrame(rows)
plan.to_csv(PLAN, index=False, encoding='utf-8-sig')

# 统一口径输入清单
with io.open(CSV_OUT, 'w', encoding='utf-8', newline='\n') as f:
    f.write('sub,lesion_path\n')
    for r in plan.itertuples():
        f.write('%s,%s\n' % (r.subject, r.mask_used))

print('总例数 %d' % len(plan))
print('复用既有 LQT 输出 : %d 例（8 月薄层批次）' % int(plan.reuse_lqt.sum()))
print('需重跑 LQT        : %d 例（9/14 重跑批次）' % int((~plan.reuse_lqt).sum()))
print('改用旧路径的例（无新掩膜）: %d 例 → %s'
      % (int((~plan.new_path_exists).sum()), list(plan[~plan.new_path_exists].subject)))
print()
sub = plan[plan.batch == '2026-08-07_thin61']
print('61 例批次：MNI/ 与 MNI_1mm_Output 掩膜一致性（用于论证复用合法）')
print('  可比对数 %d | Dice 中位 %.4f | 最小 %.4f | 完全一致（Dice>0.9999）%d 例'
      % (sub.dice_new_vs_old.notna().sum(), sub.dice_new_vs_old.median(),
         sub.dice_new_vs_old.min(), int((sub.dice_new_vs_old > 0.9999).sum())))
if (sub.dice_new_vs_old < 0.9999).any():
    print(sub[sub.dice_new_vs_old < 0.9999][['subject', 'dice_new_vs_old']].to_string(index=False))
print()
s63 = plan[plan.batch == '2026-07-09_thick63']
print('63 例批次：新掩膜 vs 旧掩膜（预期 = 旧的路线差异，中位约 0.68）')
print('  Dice 中位 %.4f | 最小 %.4f | 最大 %.4f'
      % (s63.dice_new_vs_old.median(), s63.dice_new_vs_old.min(), s63.dice_new_vs_old.max()))
print()
print('已写出:', CSV_OUT)
print('已写出:', PLAN)
