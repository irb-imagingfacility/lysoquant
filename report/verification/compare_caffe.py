"""Compare the PyTorch server output on test/TEST.tif with the segmentation made by the
original Caffe plugin (test/LQ_TEST.tif, labels 0 = background, 1 = negative, 2 = positive).

Writes compare.json and the figures compare_overview.png / compare_zoom.png / compare_margin.png
next to this script's parent folder (report/)."""
import sys, json, numpy as np, torch, tifffile
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import ndimage as ndi
R = '/Users/diegomorone/Documents/lysoquant-torch'
sys.path.insert(0, R + '/server')
from lysoquant_server.convert import load_model
from lysoquant_server import predict as P
P.log = lambda m: None
net, cfg = load_model(R + '/model/lysoquant-7-16.pt')

# same preprocessing as real.py
a = tifffile.imread(R + '/test/TEST.tif')
rgb = np.zeros((512, 512, 3), np.uint8); rgb[..., 0] = a[1]; rgb[..., 1] = a[2]
px = 1 / 8.986863; sc = px / cfg['element_size_um'][0]
fr = P.rescale_xy(P.to_frames(rgb)[0], sc, sc); fr, _, _ = P.normalize(fr, cfg['normalization_type'])
caffe = tifffile.imread(R + '/test/LQ_TEST.tif').astype(int)
assert caffe.shape == fr.shape[1:], (caffe.shape, fr.shape)

def run(tile):
    s = P.predict_frame(net, fr, tile, torch.device('cpu'))
    return s, s.argmax(0)

scores, lab = run(900)
out = dict(shape=list(caffe.shape), n_pixels=int(caffe.size))

# ---- pixel level
cm = np.array([[((caffe == i) & (lab == j)).sum() for j in range(3)] for i in range(3)])
out['confusion'] = cm.tolist(); out['diff_pixels'] = int((caffe != lab).sum()); out['accuracy'] = float(np.trace(cm) / cm.sum())
def dice(x, y): return float(2 * (x & y).sum() / (x.sum() + y.sum()))
def iou(x, y): return float((x & y).sum() / (x | y).sum())
out['per_class'] = {n: dict(dice=dice(caffe == k, lab == k), iou=iou(caffe == k, lab == k),
                            caffe_px=int((caffe == k).sum()), torch_px=int((lab == k).sum()))
                    for n, k in [('foreground', None), ('negative', 1), ('positive', 2)]
                    if k is not None} | dict(foreground=dict(dice=dice(caffe > 0, lab > 0), iou=iou(caffe > 0, lab > 0),
                                                              caffe_px=int((caffe > 0).sum()), torch_px=int((lab > 0).sum())))
# ---- where are the differing pixels?
d = caffe != lab; ys, xs = np.nonzero(d)
srt = np.sort(scores, axis=0); margin = srt[-1] - srt[-2]
edge = ndi.binary_dilation(caffe > 0) ^ ndi.binary_erosion(caffe > 0)
dist = ndi.distance_transform_edt(~edge)
clusters, ncl = ndi.label(d, structure=np.ones((3, 3)))
tm_y, tm_x = ys % 900, xs % 900
out['diff'] = dict(margin_median=float(np.median(margin[d])), margin_max=float(margin[d].max()),
                   margin_median_all=float(np.median(margin)), frac_margin_below_0p1_all=float((margin < 0.1).mean()),
                   within_2px_of_boundary=float((dist[d] <= 2).mean()), max_dist_to_boundary=float(dist[d].max()),
                   clusters=int(ncl), largest_cluster=int(np.bincount(clusters.ravel())[1:].max()),
                   within_5px_of_tile_edge=int(((np.minimum(tm_y, 900 - tm_y) <= 5) | (np.minimum(tm_x, 900 - tm_x) <= 5)).sum()),
                   within_10px_of_image_border=int(((ys < 10) | (xs < 10) | (ys > caffe.shape[0] - 11) | (xs > caffe.shape[1] - 11)).sum()),
                   frac_of_boundary_pixels=float(d.sum() / ndi.binary_dilation(edge, iterations=1).sum()))
# ---- object level (8-connected, as ImageJ's particle analyser)
obj = {}
for name, k in [('negative', 1), ('positive', 2), ('all', None)]:
    mc = (caffe > 0) if k is None else (caffe == k); mp = (lab > 0) if k is None else (lab == k)
    lc, nc = ndi.label(mc, structure=np.ones((3, 3))); lp, npp = ndi.label(mp, structure=np.ones((3, 3)))
    exact = 0; ious = []
    for i, sl in enumerate(ndi.find_objects(lc), 1):
        m = lc[sl] == i; ov = lp[sl][m]; ov = ov[ov > 0]
        if ov.size == 0: ious.append(0.); continue
        j = np.bincount(ov).argmax(); mj = lp == j
        ious.append(float(((lc == i) & mj).sum() / ((lc == i) | mj).sum()))
        exact += int(ious[-1] == 1.0)
    ac, ap = np.bincount(lc.ravel())[1:], np.bincount(lp.ravel())[1:]
    obj[name] = dict(caffe=int(nc), torch=int(npp), identical=exact, mean_iou=float(np.mean(ious)), min_iou=float(np.min(ious)),
                     area_caffe=int(ac.sum()), area_torch=int(ap.sum()))
out['objects'] = obj
# ---- how much does the tiling matter?
sweep = []
for T in (388, 516, 900, 1156):
    l = lab if T == 900 else run(T)[1]
    sweep.append(dict(tile=T, tiles=int(np.prod([P.tile_layout(n, T, net)[1] for n in caffe.shape])),
                      diff_vs_caffe=int((l != caffe).sum()), frac=float((l != caffe).mean())))
    print(sweep[-1], flush=True)
out['tile_sweep'] = sweep
json.dump(out, open(R + '/report/verification/compare.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'confusion'}, indent=1))

# ---- figures
NEG, POS = '#ffc000', '#ff00ff'
cmap = ListedColormap(['black', NEG, POS])
fig, ax = plt.subplots(1, 3, figsize=(13, 4.5))
ax[0].imshow(caffe, cmap=cmap, vmin=0, vmax=2, interpolation='nearest'); ax[0].set_title('Caffe (original plugin)')
ax[1].imshow(lab, cmap=cmap, vmin=0, vmax=2, interpolation='nearest'); ax[1].set_title('PyTorch (this server)')
ax[2].imshow(caffe > 0, cmap='gray_r', alpha=.25); ax[2].scatter(xs, ys, s=22, facecolors='none', edgecolors='red', linewidths=.7)
ax[2].set_title(f'{d.sum()} differing pixels of {d.size:,}'); ax[2].set_xlim(0, caffe.shape[1]); ax[2].set_ylim(caffe.shape[0], 0)
for x in ax: x.axis('off')
fig.tight_layout(); fig.savefig(R + '/report/compare_overview.png', dpi=130, bbox_inches='tight'); plt.close(fig)

# zoom: the two 100x100 windows holding most differing pixels
H = 50; cnt = ndi.uniform_filter(d.astype(float), 2 * H, mode='constant') * (2 * H) ** 2
wins = []
c2 = cnt.copy()
for _ in range(2):
    y, x = np.unravel_index(c2.argmax(), c2.shape); y = int(np.clip(y, H, caffe.shape[0] - H)); x = int(np.clip(x, H, caffe.shape[1] - H))
    wins.append((y, x)); c2[max(y - 2 * H, 0):y + 2 * H, max(x - 2 * H, 0):x + 2 * H] = 0
fig, ax = plt.subplots(2, 3, figsize=(10, 6.7))
for r, (y, x) in enumerate(wins):
    sl = (slice(y - H, y + H), slice(x - H, x + H)); ext = (x - H, x + H, y + H, y - H)
    ax[r, 0].imshow(caffe[sl], cmap=cmap, vmin=0, vmax=2, interpolation='nearest', extent=ext)
    ax[r, 1].imshow(lab[sl], cmap=cmap, vmin=0, vmax=2, interpolation='nearest', extent=ext)
    ax[r, 2].imshow(caffe[sl] > 0, cmap='gray_r', alpha=.3, interpolation='nearest', extent=ext)
    yy, xx = np.nonzero(d[sl]); ax[r, 2].scatter(xx + x - H + .5, yy + y - H + .5, s=40, facecolors='none', edgecolors='red', linewidths=.8)
    ax[r, 2].set_xlim(ext[0], ext[1]); ax[r, 2].set_ylim(ext[2], ext[3])
    ax[r, 2].set_title(f'{len(yy)} differing px'); ax[r, 0].set_title('Caffe'); ax[r, 1].set_title('PyTorch')
    for q in ax[r]: q.set_xticks([]); q.set_yticks([])
fig.tight_layout(); fig.savefig(R + '/report/compare_zoom.png', dpi=130, bbox_inches='tight'); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.5, 3.4))
bins = np.logspace(-4, 1.6, 60)
ax.hist(np.clip(margin[~d], 1e-4, None).ravel(), bins=bins, weights=np.full((~d).sum(), 1 / (~d).sum()), color='0.6', label='agreeing pixels')
ax.hist(margin[d], bins=bins, weights=np.full(d.sum(), 1 / d.sum()), color='red', alpha=.7, label='differing pixels')
ax.set_xscale('log'); ax.set_xlabel('score margin (top-1 minus top-2 logit)'); ax.set_ylabel('fraction of pixels'); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(R + '/report/compare_margin.png', dpi=130); plt.close(fig)
