import sys, json, time
import numpy as np, h5py, torch
sys.path.insert(0, '/Users/diegomorone/Documents/lysoquant-torch/server')
from ptparse import load
from lysoquant_server.convert import load_model
R = '/Users/diegomorone/Documents/lysoquant-torch/model/'
res = {}

# ---- test-phase layer list from the prototxt
proto = load('model.prototxt')
layers = [L for L in proto['layer'] if not L.get('include') or L['include'][0]['phase'][0] != 'TRAIN']
res['n_layers_test'] = len(layers)

net, cfg = load_model(R + 'lysoquant-7-16.pt')
res['config'] = {k: v for k, v in cfg.items()}

# ---- 1. structure: prototxt vs module graph
conv_like = [L for L in layers if L['type'][0] in ('Convolution', 'Deconvolution')]
mods = []
for b in net.down: mods += [b[0], b[2]]
# execution order in forward(): down convs, then per level: upconv, up.0, up.2
for uc, ub in zip(net.upconv, net.up): mods += [uc, ub[0], ub[2]]
mods.append(net.score)
assert len(mods) == len(conv_like) == 38, (len(mods), len(conv_like))
rows = []
ok = True
for L, m in zip(conv_like, mods):
    p = L['convolution_param'][0]
    k = [int(x) for x in p['kernel_size']]; k = k * (2 // len(k)) if len(k) == 1 else k
    st = int(p.get('stride', ['1'])[0]); pad = int(p.get('pad', ['0'])[0])
    for key in ('group', 'dilation', 'bias_term'):
        assert key not in p, key
    isdec = L['type'][0] == 'Deconvolution'
    good = (isdec == isinstance(m, torch.nn.ConvTranspose2d) and tuple(m.kernel_size) == tuple(k[:2] if len(k)>1 else k*2)
            and m.stride[0] == st and m.padding[0] == pad
            and (m.weight.shape[1] if isdec else m.weight.shape[0]) == int(p['num_output'][0]))
    ok &= good
    rows.append((L['name'][0], L['type'][0], str(k), st, pad, int(p['num_output'][0]), good))
res['structure_ok'] = bool(ok); res['structure_rows'] = rows
slopes = {s['negative_slope'][0] for L in layers if L['type'][0]=='ReLU' for s in L['relu_param']}
res['slopes'] = sorted(slopes)
res['n_relu'] = sum(L['type'][0]=='ReLU' for L in layers)

# ---- 2. weights: bit-exact, mapping by execution order (independent of convert.py's name map)
f = h5py.File(R + '2D-7-16.caffemodel.h5')
maxdiff = 0; nparam = 0; exact = True
for L, m in zip(conv_like, mods):
    w = np.asarray(f['data'][L['name'][0]]['0'], np.float32); b = np.asarray(f['data'][L['name'][0]]['1'], np.float32)
    exact &= np.array_equal(w, m.weight.detach().numpy()) and np.array_equal(b, m.bias.detach().numpy())
    nparam += w.size + b.size
res['weights_bit_exact'] = bool(exact); res['n_params'] = int(nparam)
res['h5_param_layers'] = len(f['data'])
res['pt_params'] = sum(p.numel() for p in net.parameters())

# ---- 3. independent NumPy interpreter of the prototxt (float64)
from numpy.lib.stride_tricks import sliding_window_view as swv
def conv(x, w, b):
    kh, kw = w.shape[2:]; H, W = x.shape[1] - kh + 1, x.shape[2] - kw + 1
    out = np.zeros((w.shape[0], H, W))
    for i in range(kh):
        for j in range(kw):
            out += np.tensordot(w[:, :, i, j], x[:, i:i + H, j:j + W], axes=([1], [0]))
    return out + b[:, None, None]
def deconv(x, w, b, s):                              # w: (Cin, Cout, k, k), stride s == k
    C, H, W = x.shape; k = w.shape[2]
    out = np.zeros((w.shape[1], H * s, W * s))
    for i in range(k):
        for j in range(k):
            out[:, i::s, j::s] += np.tensordot(w[:, :, i, j], x, axes=([0], [0]))
    return out + b[:, None, None]
def pool(x):
    C, H, W = x.shape
    return x.reshape(C, H // 2, 2, W // 2, 2).max(axis=(2, 4))
def crop(x, h, w):
    dy, dx = (x.shape[1] - h) // 2, (x.shape[2] - w) // 2
    return x[:, dy:dy + h, dx:dx + w]
def run_proto(x):
    blobs = {'data3': x.astype(np.float64)}
    for L in layers:
        t, name = L['type'][0], L['name'][0]
        bot = [blobs[b] for b in L.get('bottom', [])]; top = L['top'][0]
        if t == 'ValueTransformation':
            o = bot[0] + float(L['value_transformation_param'][0]['offset'][0]['v'][0])
        elif t == 'Convolution':
            o = conv(bot[0], np.asarray(f['data'][name]['0'], np.float64), np.asarray(f['data'][name]['1'], np.float64))
        elif t == 'Deconvolution':
            st = int(L['convolution_param'][0]['stride'][0])
            o = deconv(bot[0], np.asarray(f['data'][name]['0'], np.float64), np.asarray(f['data'][name]['1'], np.float64), st)
        elif t == 'ReLU':
            s = float(L['relu_param'][0]['negative_slope'][0]); o = np.where(bot[0] > 0, bot[0], s * bot[0])
        elif t == 'Pooling':
            o = pool(bot[0])
        elif t == 'Concat':
            a, c = bot; o = np.concatenate([a, crop(c, a.shape[1], a.shape[2])], axis=0)  # caffe_unet Concat crops to the smallest blob
        else:
            raise ValueError(t)
        blobs[top] = o
    return blobs['score']

out_sizes = {n: net.output_shape(n) for n in range(1, 2000)}
valid = [n for n, o in out_sizes.items() if o]
res['min_input'], res['min_output'] = valid[0], out_sizes[valid[0]]
res['io_offset'] = valid[0] - out_sizes[valid[0]]
res['valid_step'] = valid[1] - valid[0]
res['proto_pad'] = {k: f2[()].tolist() for f2, k in [(h5py.File(R+'2D-7-16.modeldef.h5')['unet_param'][k], k) for k in ('padInput','padOutput','downsampleFactor')]}

rng = np.random.default_rng(0)
tests = []
for n in (valid[0], valid[1]):
    if n==valid[1] and False: continue
    for kind in ('uniform', 'blobs'):
        if kind == 'uniform':
            x = rng.random((3, n, n))
        else:   # smooth blobs + noise, more realistic dynamic range
            yy, xx = np.mgrid[:n, :n]; x = np.stack([np.exp(-(((yy-rng.integers(n))**2+(xx-rng.integers(n))**2)/(2*(n/8)**2))) for _ in range(3)])
            x = np.clip(x + 0.05 * rng.standard_normal(x.shape), 0, 1)
        t0 = time.time(); ref = run_proto(x); t1 = time.time()
        with torch.no_grad():
            torch.backends.cudnn.allow_tf32 = False
            got = net(torch.from_numpy(x.astype(np.float32))[None])[0].numpy()
        # float64 torch to separate precision from structure
        net64 = net.double()
        with torch.no_grad(): got64 = net64(torch.from_numpy(x)[None])[0].numpy()
        net.float()
        amax = np.abs(ref).max()
        tests.append(dict(n=n, kind=kind, out=list(ref.shape), t_ref=t1 - t0,
                          max_abs32=float(np.abs(ref - got).max()), max_abs64=float(np.abs(ref - got64).max()),
                          scale=float(amax), rel32=float(np.abs(ref - got).max() / amax),
                          argmax_agree=float((ref.argmax(0) == got.argmax(0)).mean()),
                          classes=np.bincount(ref.argmax(0).ravel(), minlength=3).tolist()))
        print(tests[-1], flush=True)
res['tests'] = tests
json.dump(res, open('res.json', 'w'), indent=1, default=str)
print(json.dumps({k: v for k, v in res.items() if k not in ('structure_rows', 'tests')}, indent=1, default=str))
