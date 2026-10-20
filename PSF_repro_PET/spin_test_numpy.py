"""Pure-numpy spin permutation test for cortical surface maps.

Drop-in alternative to enigmatoolbox.permutation_testing.spin_test that
avoids importing vtk (load_fsa5), which crashes on many Windows installs.
Reads fsa5 sphere .gii and aparc .annot with nibabel instead.

The three rotation functions (centroid_extraction_sphere, rotate_parcellation,
perm_sphere_p) are verbatim from enigmatoolbox (MICA-MNI, @saratheriver),
released for reuse; see Alexander-Bloch et al. 2018, NeuroImage 178:540-51 and
Vasa et al. 2017, Cereb Cortex 28(1):281-294.

Usage:
    from spin_test_numpy import spin_test_numpy
    p_spin, r_dist = spin_test_numpy(map1, map2, surface_dir, annot_dir, n_rot=1000)
"""
import os
import warnings
import numpy as np
import pandas as pd
import nibabel as nb
from scipy.spatial.distance import cdist


def load_fsa5_sphere_numpy(surface_dir):
    """Read fsa5 sphere coordinates via nibabel (no vtk).
    Returns (sphere_lh, sphere_rh), each (n_vertices, 3)."""
    def _read(hemi):
        g = nb.load(os.path.join(surface_dir, f'fsa5_sphere_{hemi}.gii'))
        return g.darrays[0].data
    return _read('lh'), _read('rh')


# ---- verbatim from enigmatoolbox.permutation_testing ----
def centroid_extraction_sphere(sphere_coords, annotfile, ventricles=False):
    if "aparc_aseg" not in annotfile:
        labels, ctab, names = nb.freesurfer.io.read_annot(annotfile, orig_ids=True)
        centroid = np.empty((0, 3))
        for ic in range(ctab.shape[0]):
            if not names[ic].decode("utf-8") == 'unknown' and not names[ic].decode("utf-8") == 'corpus'\
                and not names[ic].decode("utf-8") == 'corpuscallosum' and not names[ic].decode("utf-8") == 'medialwall'\
                and not names[ic].decode("utf-8") == 'Background+FreeSurfer_Defined_Medial_Wall':
                label = ctab[ic, -1]
                centroid = np.vstack((centroid, np.array(np.mean(sphere_coords[labels == label, :], axis=0))))
    elif "aparc_aseg" in annotfile and ventricles is not True:
        annot = pd.read_csv(annotfile)
        centroid = np.empty((0, 3))
        for ic in range(annot['structure'][:44].shape[0]):
            if not annot['structure'][:44][ic] == "'unknown'" and not annot['structure'][:44][ic] == "'corpuscallosum'" and not annot['structure'][:44][ic] == "'vent'":
                label = annot['label'][ic]
                centroid = np.vstack((centroid, np.array(np.mean(sphere_coords[annot[annot['label_annot'] ==
                                                                                     label].index.values, :], axis=0))))
    elif "aparc_aseg" in annotfile and ventricles is True:
        annot = pd.read_csv(annotfile)
        centroid = np.empty((0, 3))
        for ic in range(annot['structure'][:44].shape[0]):
            if not annot['structure'][:44][ic] == "'unknown'" and not annot['structure'][:44][ic] == "'corpuscallosum'":
                label = annot['label'][ic]
                centroid = np.vstack((centroid, np.array(np.mean(sphere_coords[annot[annot['label_annot'] ==
                                                                                     label].index.values, :], axis=0))))
    return centroid


def rotate_parcellation(coord_l, coord_r, nrot=1000, seed=0):
    """Generate spin-rotation permutation index matrix.

    NOTE: a fixed ``seed`` (default=0) is used so the rotation matrix is
    DETERMINISTIC across runs/processes. Previously this used the global
    ``np.random`` state with no seed, which made every downstream spin p-value
    non-reproducible. Pass ``seed=None`` to restore the old unseeded behaviour.
    """
    warnings.filterwarnings('ignore')
    if coord_l.shape[1] != 3 or coord_r.shape[1] != 3:
        coord_l = np.transpose(coord_l)
        coord_r = np.transpose(coord_r)
    nroi_l = coord_l.shape[0]
    nroi_r = coord_r.shape[0]
    nroi = nroi_l + nroi_r
    perm_id = np.zeros((nroi, nrot))
    r = 0
    c = 0
    I1 = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
    _rng = np.random if seed is None else np.random.default_rng(seed)
    while (r < nrot):
        A = _rng.normal(loc=0, scale=1, size=(3, 3))
        TL, temp = np.linalg.qr(A)
        TL = np.matmul(TL, np.diag(np.sign(np.diag(temp))))
        if np.linalg.det(TL) < 0:
            TL[:, 0] = -TL[:, 0]
        TR = np.matmul(np.matmul(I1, TL), I1)
        coord_l_rot = np.matmul(coord_l, TL)
        coord_r_rot = np.matmul(coord_r, TR)
        dist_l = cdist(coord_l, coord_l_rot)
        dist_r = cdist(coord_r, coord_r_rot)
        temp_dist_l = dist_l
        rot_l = []
        ref_l = []
        for _ in range(nroi_l):
            ref_ix = np.argwhere(np.nanmin(temp_dist_l, axis=1) == np.nanmax(np.nanmin(temp_dist_l, axis=1)))[0]
            rot_ix = np.argwhere(temp_dist_l[ref_ix, :][0] == np.nanmin(temp_dist_l[ref_ix, :], axis=1))[0]
            ref_l = np.append(ref_l, ref_ix)
            rot_l = np.append(rot_l, rot_ix)
            temp_dist_l[:, rot_ix] = np.nan
            temp_dist_l[ref_ix, :] = np.nan
        temp_dist_r = dist_r
        rot_r = []
        ref_r = []
        for _ in range(nroi_r):
            ref_ix = np.argwhere(np.nanmin(temp_dist_r, axis=1) == np.nanmax(np.nanmin(temp_dist_r, axis=1)))[0]
            rot_ix = np.argwhere(temp_dist_r[ref_ix, :][0] == np.nanmin(temp_dist_r[ref_ix, :], axis=1))[0]
            ref_r = np.append(ref_r, ref_ix)
            rot_r = np.append(rot_r, rot_ix)
            temp_dist_r[:, rot_ix] = np.nan
            temp_dist_r[ref_ix, :] = np.nan
        ref_lr = np.append(ref_l, nroi_l + ref_r)
        rot_lr = np.append(rot_l, nroi_l + rot_r)
        ix = np.argsort(ref_lr)
        rot_lr_sort = rot_lr[ix]
        if np.all(rot_lr_sort == range(nroi)) is not True:
            perm_id[:, r] = rot_lr_sort
            r = r + 1
        elif np.all(rot_lr_sort == range(nroi)) is True:
            c = c + 1
    return perm_id


def perm_sphere_p(x, y, perm_id, corr_type='pearson', null_dist=False):
    nroi = perm_id.shape[0]
    nperm = perm_id.shape[1]
    rho_emp = pd.Series(x).corr(pd.Series(y), method=corr_type)
    x_perm = np.empty((0, nroi))
    y_perm = np.empty((0, nroi))
    for rr in range(nperm):
        x_perm2 = []
        y_perm2 = []
        for ii in range(nroi):
            x_perm2 = np.append(x_perm2, x[int(perm_id[ii, rr])])
            y_perm2 = np.append(y_perm2, y[int(perm_id[ii, rr])])
        x_perm = np.vstack((x_perm, x_perm2))
        y_perm = np.vstack((y_perm, y_perm2))
    x_perm = np.transpose(x_perm)
    y_perm = np.transpose(y_perm)
    rho_null_xy = []
    rho_null_yx = []
    for rr in range(nperm):
        rho_null_xy = np.append(rho_null_xy, pd.Series(x_perm[:, rr]).corr(pd.Series(y), method=corr_type))
        rho_null_yx = np.append(rho_null_yx, pd.Series(x).corr(pd.Series(y_perm[:, rr]), method=corr_type))
    if rho_emp >= 0:
        p_perm_xy = np.sum((rho_null_xy > rho_emp).astype(int)) / nperm
        p_perm_yx = np.sum((rho_null_yx > rho_emp).astype(int)) / nperm
    elif rho_emp < 0:
        p_perm_xy = np.sum((rho_null_xy < rho_emp).astype(int)) / nperm
        p_perm_yx = np.sum((rho_null_yx < rho_emp).astype(int)) / nperm
    p_perm = (p_perm_xy + p_perm_yx) / 2
    if null_dist is True:
        return p_perm, np.append(rho_null_xy, rho_null_yx)
    elif null_dist is not True:
        return p_perm


def spin_test_numpy(map1, map2, surface_dir, annot_dir,
                    parcellation_name='aparc', n_rot=1000,
                    corr_type='pearson'):
    """Spin permutation test, vtk-free. Mirrors enigmatoolbox.spin_test.

    map1, map2 : 1D arrays (length 68 for aparc cortical), same region order.
    surface_dir : dir containing fsa5_sphere_lh.gii / fsa5_sphere_rh.gii
                  (enigmatoolbox/datasets/surfaces).
    annot_dir   : dir containing fsa5_lh_aparc.annot / fsa5_rh_aparc.annot
                  (enigmatoolbox/permutation_testing/annot).
    Returns (p_spin, r_dist) where r_dist is the null distribution (n_rot*2,).
    """
    sphere_lh, sphere_rh = load_fsa5_sphere_numpy(surface_dir)
    annotfile_lh = os.path.join(annot_dir, f'fsa5_lh_{parcellation_name}.annot')
    annotfile_rh = os.path.join(annot_dir, f'fsa5_rh_{parcellation_name}.annot')
    lh_centroid = centroid_extraction_sphere(sphere_lh, annotfile_lh)
    rh_centroid = centroid_extraction_sphere(sphere_rh, annotfile_rh)
    perm_id = rotate_parcellation(lh_centroid, rh_centroid, nrot=n_rot)
    if isinstance(map1, (pd.DataFrame, pd.Series)):
        map1 = map1.to_numpy()
    if isinstance(map2, (pd.DataFrame, pd.Series)):
        map2 = map2.to_numpy()
    return perm_sphere_p(map1, map2, perm_id, corr_type, null_dist=True)


if __name__ == "__main__":
    # smoke test needs a real enigmatoolbox install for surface_dir/annot_dir.
    # Replace the two paths below with your local enigmatoolbox package dirs,
    # then run:  python spin_test_numpy.py
    import sys
    if len(sys.argv) >= 3:
        sdir, adir = sys.argv[1], sys.argv[2]
        rng = np.random.default_rng(0)
        a = rng.normal(0, 1, 68)
        b = rng.normal(0, 1, 68)
        p, rd = spin_test_numpy(a, b, sdir, adir, n_rot=200)
        print(f"smoke test p={p:.3f} null_len={len(rd)}")
    else:
        print("usage: python spin_test_numpy.py <surface_dir> <annot_dir>")
