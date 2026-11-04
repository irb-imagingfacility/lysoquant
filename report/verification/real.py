import sys, json, time, numpy as np, torch, tifffile
sys.path.insert(0,'/Users/diegomorone/Documents/lysoquant-torch/server')
from lysoquant_server.convert import load_model
from lysoquant_server import predict as P
from refnet import run_proto
P.log = lambda m: None
net, cfg = load_model('/Users/diegomorone/Documents/lysoquant-torch/model/lysoquant-7-16.pt')
a = tifffile.imread('/Users/diegomorone/Documents/lysoquant-torch/test/TEST.tif')
rgb = np.zeros((512,512,3), np.uint8); rgb[...,0]=a[1]; rgb[...,1]=a[2]  # ImageJ ch2=protein (R), ch3=lysosome (G)
px = 1/8.986863; el = [px, px]; sc = px/cfg['element_size_um'][0]
frames = P.to_frames(rgb)
fr = P.rescale_xy(frames[0], sc, sc); fr, scale, off = P.normalize(fr, cfg['normalization_type'])
H, W = fr.shape[1:]; print('rescaled', fr.shape, scale, off, flush=True)
# full pipeline
t0=time.time(); scores = P.predict_frame(net, fr, 900, torch.device('cpu')); lab = scores.argmax(0).astype(np.uint16)
print('pipeline s', time.time()-t0, np.bincount(lab.ravel(), minlength=3), flush=True)
np.save('lab.npy', lab); np.save('rgb.npy', rgb)
# tile-level reference check on real, preprocessed, mirror-padded data
pad = 764; T = 132
rows = P.mirror_indices(-pad, H+pad, H); cols = P.mirror_indices(-pad, W+pad, W)
canvas = fr[:, rows[:,None], cols]
res = []
for (oy, ox) in [(1000, 1000), (300, 1400), (1900, 500), (0, 0)]:
    x = canvas[:, oy:oy+T+2*pad, ox:ox+T+2*pad]
    ref = run_proto(x)
    with torch.no_grad():
        g32 = net(torch.from_numpy(np.ascontiguousarray(x))[None])[0].numpy()
        g64 = net.double()(torch.from_numpy(np.ascontiguousarray(x, dtype=np.float64))[None])[0].numpy(); net.float()
    r = dict(origin=[oy,ox], max_abs32=float(np.abs(ref-g32).max()), max_abs64=float(np.abs(ref-g64).max()),
             scale=float(np.abs(ref).max()), argmax_agree=float((ref.argmax(0)==g32.argmax(0)).mean()),
             classes=np.bincount(ref.argmax(0).ravel(), minlength=3).tolist(),
             classes_torch=np.bincount(g32.argmax(0).ravel(), minlength=3).tolist())
    print(r, flush=True); res.append(r)
json.dump(dict(shape=[H,W], tiles=res, classes=np.bincount(lab.ravel(), minlength=3).tolist()), open('real.json','w'))
