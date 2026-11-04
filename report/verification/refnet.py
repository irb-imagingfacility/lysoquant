import numpy as np, h5py
from ptparse import load
R='/Users/diegomorone/Documents/lysoquant-torch/model/'
proto=load('model.prototxt')
layers=[L for L in proto['layer'] if not L.get('include') or L['include'][0]['phase'][0]!='TRAIN']
f=h5py.File(R+'2D-7-16.caffemodel.h5')
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

