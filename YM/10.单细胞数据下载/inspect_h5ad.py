# 探查 CELLxGENE h5ad 结构
import anndata as ad
import numpy as np

a = ad.read_h5ad("periodontal_atlas.h5ad")
print("==== 1. 维度 (细胞 x 基因) ====")
print(a.shape)

print("\n==== 2. X 矩阵 ====")
print("type:", type(a.X), "dtype:", a.X.dtype)
X = a.X
sample = np.asarray(X[:5000].toarray() if hasattr(X, "toarray") else X[:5000])
print("前5000细胞  min/median/max: %.3f / %.3f / %.3f" %
      (sample.min(), np.median(sample), sample.max()))
print("是否已归一化(最大值<50):", sample.max() < 50)

print("\n==== 3. obs 列 ====")
print(list(a.obs.columns))
for c in a.obs.columns:
    s = a.obs[c]
    nu = s.nunique()
    print("\n-- %s  (%s)  n_unique=%d" % (c, s.dtype, nu))
    if nu <= 40:
        print(s.value_counts().to_string())
    else:
        print("  取值过多, 前10:", list(s.dropna().unique()[:10]))

print("\n==== 4. var（基因）====")
print("var 列:", list(a.var.columns))
print("前10基因:", list(a.var_names[:10]))

print("\n==== 5. obsm（降维）====")
print(list(a.obsm_keys()))
for k in a.obsm_keys():
    print(" ", k, a.obsm[k].shape)

print("\n==== 6. raw / layers ====")
print("raw 存在:", a.raw is not None)
if a.raw is not None:
    print("  raw shape:", a.raw.shape)
print("layers:", list(a.layers.keys()))

print("\n==== 7. uns（前30键）====")
print(list(a.uns.keys())[:30])
