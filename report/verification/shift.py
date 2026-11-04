import sys, json, numpy as np, torch
sys.path.insert(0,'/Users/diegomorone/Documents/lysoquant-torch/server')
from lysoquant_server.convert import load_model
net, cfg = load_model('/Users/diegomorone/Documents/lysoquant-torch/model/lysoquant-7-16.pt')
rng = np.random.default_rng(2)
N = 1788; S = 260
# smooth-ish random RGB canvas with bright blobs
yy, xx = np.mgrid[:N+256, :N+256]; canvas = np.zeros((3, N+256, N+256), np.float32)
for _ in range(120):
    cy, cx, r = rng.integers(0, N+256), rng.integers(0, N+256), rng.integers(6, 30)
    canvas[rng.integers(3)] += np.exp(-(np.hypot(yy-cy, xx-cx)/r)**2)
canvas = (canvas - canvas.min())/(canvas.max()-canvas.min())
def run(dx):
    x = torch.from_numpy(np.ascontiguousarray(canvas[:, :N, dx:dx+N]))[None]
    with torch.no_grad(): return net(x)[0].numpy()
base = run(0); out = {}
for dx in (128, 256, 4, 64, 132):
    o = run(dx)
    # output col j of crop dx corresponds to output col j+dx of crop 0
    a = base[:, :, dx:]; b = o[:, :, :S-dx]
    out[dx] = dict(max_abs=float(np.abs(a-b).max()), argmax_agree=float((a.argmax(0)==b.argmax(0)).mean()), scale=float(np.abs(base).max()))
    print(dx, out[dx], flush=True)
json.dump(out, open('shift.json','w'))
