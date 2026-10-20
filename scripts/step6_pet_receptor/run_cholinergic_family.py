# -*- coding: utf-8 -*-
"""R3：胆碱能家族内部一致性检查（三张胆碱能图谱）
复现受体管线，检查 flubatine(α4β2) / feobv(VAChT) / lsn3172176(M1) 三者之间的
空间相关，以及彼此控制后的偏相关与 spin p。"""
import os, sys
import numpy as np, pandas as pd, nibabel as nb
from scipy.ndimage import affine_transform
BASE='F:/dcmtonii/run62_crosssectional'; PET=BASE+'/PSF_repro_PET'
SAL=BASE+'/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
ATLAS=BASE+'/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.nii.gz'
NM=PET+'/data/neuromaps_annotations'; SEED=0; NROT=1000
sal=pd.read_csv(SAL); lv1=sal['salience'].values.astype(float); parcels=sal['region'].tolist()
ai=nb.load(ATLAS); ad=ai.get_fdata(); aa=ai.affine
def parc(path):
    img=nb.load(path); rec=img.get_fdata().squeeze()
    if rec.shape!=ad.shape or not np.allclose(img.affine,aa):
        c=np.linalg.inv(aa)@img.affine
        ar=affine_transform(ad,c[:3,:3],offset=c[:3,3],output_shape=rec.shape,order=0).astype(int)
    else: ar=ad.astype(int)
    v=np.full(100,np.nan)
    for lab in range(1,101):
        m=(ar==lab)&~np.isnan(rec)
        if m.sum()>0: v[lab-1]=np.nanmean(rec[m])
    return np.where(np.isnan(v),np.nanmean(v),v)
CH={'flubatine (α4β2)':('hillmer2016','flubatine','1mm'),
    'feobv (VAChT)':('aghourian2017','feobv','1mm'),
    'lsn3172176 (M1)':('naganawa2020','lsn3172176','1mm')}
M={k:parc(os.path.join(NM,s,t,'MNI152',f'source-{s}_desc-{t}_space-MNI152_res-{r}_feature.nii.gz'))
   for k,(s,t,r) in CH.items()}
sys.path.insert(0, PET)
from spin_test_numpy import rotate_parcellation, perm_sphere_p
sd=f'{PET}/assets'
sl=nb.load(f'{sd}/fsa5_sphere_lh.gii').darrays[0].data; sr=nb.load(f'{sd}/fsa5_sphere_rh.gii').darrays[0].data
def coords_by_order(annot_f, sphere):
    # 与 run_receptor_discon_cortex100_residY.py 逐字一致
    labels, ctab, names = nb.freesurfer.io.read_annot(annot_f)
    names = [nm.decode() for nm in names]
    out = {}
    for i, nm in enumerate(names):
        if nm in ('unknown', 'medialwall'):
            continue
        m = labels == i
        if m.sum() == 0:
            continue
        out[nm.replace('7Networks_', '')] = sphere[m].mean(0)
    return out
Cl=coords_by_order(f'{sd}/schaefer100_fsa5/lh.Schaefer2018_100Parcels_7Networks_order.annot',sl)
Cr=coords_by_order(f'{sd}/schaefer100_fsa5/rh.Schaefer2018_100Parcels_7Networks_order.annot',sr)
cl=[];cr=[]
for p in parcels:
    hemi=p.split('_')[0]
    (cl if hemi=='LH' else cr).append((Cl if hemi=='LH' else Cr)[p])
pid=rotate_parcellation(np.array(cl),np.array(cr),nrot=NROT,seed=SEED)
def rot(v,rr): return v[pid[:,rr].astype(int)]
def pc(x,y,Z):
    Zc=np.column_stack([np.ones(len(x)),Z]); 
    def rs(v):
        b,_,_,_=np.linalg.lstsq(Zc,v,rcond=None); return v-Zc@b
    return float(np.corrcoef(rs(x),rs(y))[0,1])
def pincl(v,obs):
    # 含入式单侧计数 (1+k)/(N+1)（与论文数值口径一致；官方严格版为 k/N）
    if obs>=0: return (1+int((np.asarray(v)>=obs).sum()))/(NROT+1)
    return (1+int((np.asarray(v)<=obs).sum()))/(NROT+1)
def spinp(x,y,Z,obs):
    a=np.array([pc(rot(x,r),y,Z) for r in range(NROT)]); b=np.array([pc(x,rot(y,r),Z) for r in range(NROT)])
    return float((pincl(a,obs)+pincl(b,obs))/2)
print('=== 三张胆碱能图谱两两空间相关（100 皮层区）===')
ks=list(CH)
for i in range(len(ks)):
    for j in range(i+1,len(ks)):
        print(f'  {ks[i]:18s} vs {ks[j]:18s} r = {np.corrcoef(M[ks[i]],M[ks[j]])[0,1]:+.3f}')
print()
print('=== 与 LV1 loading 的关联 ===')
print('| 图谱 | 简单 r | spin p | 控另两张胆碱能图谱的偏 r | spin p |')
print('|---|---|---|---|---|')
for k in ks:
    x=M[k]; Z=np.column_stack([M[o] for o in ks if o!=k])
    rs=float(np.corrcoef(lv1,x)[0,1]); _,nsp=perm_sphere_p(lv1,x,pid,corr_type='pearson',null_dist=True)
    ps=float((pincl(nsp[:NROT],rs)+pincl(nsp[NROT:],rs))/2)
    ob=pc(x,lv1,Z); pb=spinp(x,lv1,Z,ob)
    print(f'| {k} | {rs:+.3f} | {ps:.4f} | {ob:+.3f} | {pb:.4f} |')
