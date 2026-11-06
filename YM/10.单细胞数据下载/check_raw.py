import anndata as ad
import numpy as np
a = ad.read_h5ad("periodontal_atlas.h5ad")
print("X (归一化) 前2000细胞 max: %.2f" % np.asarray(a.X[:2000].toarray()).max())
if a.raw is not None:
    rx = a.raw.X
    s = np.asarray(rx[:2000].toarray() if hasattr(rx, "toarray") else rx[:2000])
    print("raw.X min/max: %.2f / %.2f" % (s.min(), s.max()))
    print("raw.X 是否整数: %s" % np.allclose(s, np.round(s)))
    print("raw.var 列:", list(a.raw.var.columns))
print("var feature_name 示例:", list(a.var['feature_name'][:8]))
