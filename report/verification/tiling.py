import sys, json, numpy as np, torch
sys.path.insert(0,'/Users/diegomorone/Documents/lysoquant-torch/server')
from lysoquant_server.convert import load_model
from lysoquant_server import predict as P
P.log = lambda m: None
net, cfg = load_model('/Users/diegomorone/Documents/lysoquant-torch/model/lysoquant-7-16.pt')
rng = np.random.default_rng(1)
H, W = 200, 230
yy, xx = np.mgrid[:H, :W]; img = np.zeros((3, H, W), np.float32)
for _ in range(14):                         # ring-like blobs in G, filled blobs in R/B
    cy, cx, r = rng.integers(20, H-20), rng.integers(20, W-20), rng.integers(6, 16)
    d = np.hypot(yy-cy, xx-cx)
    img[1] += np.exp(-((d-r)/2.0)**2); img[0] += (d < r*0.6); img[2] += 0.5*(d < r)
img += 0.03*rng.standard_normal(img.shape).astype(np.float32)
img, _, _ = P.normalize(img, cfg['normalization_type'])
dev = torch.device('cpu')
one  = P.predict_frame(net, img, 900, dev)       # single tile
tiled = P.predict_frame(net, img, 132, dev)      # 2x2 tiles
d = np.abs(one - tiled)
r = dict(shape=list(one.shape), max_abs=float(d.max()), scale=float(np.abs(one).max()),
         argmax_agree=float((one.argmax(0) == tiled.argmax(0)).mean()),
         classes=np.bincount(one.argmax(0).ravel(), minlength=3).tolist(),
         ntiles=[P.tile_layout(H,132,net), P.tile_layout(W,132,net)])
print(r); json.dump(r, open('tiling.json','w'))
