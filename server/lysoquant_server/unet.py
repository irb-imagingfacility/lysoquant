"""
PyTorch re-implementation of the 2D U-Net used by caffe_unet
(Falk et al., Nature Methods 2019) as trained for LysoQuant.

The layer layout mirrors the Caffe prototxt one-to-one so that the
converted Caffe weights can be loaded without any reshaping:

  * valid (unpadded) 3x3 convolutions followed by LeakyReLU(0.1)
  * 2x2 max pooling, stride 2
  * 2x2 up-convolutions (stride 2) followed by LeakyReLU(0.1)
  * concatenation [upsampled, center-cropped skip connection]
  * final 1x1 convolution to the class scores
  * the input is shifted by -0.5 (Caffe ValueTransformation layer)
"""

import itertools

import torch
import torch.nn as nn
import torch.nn.functional as F


def _conv_block(in_ch, out_ch, slope):
    return nn.Sequential(
        nn.Conv2d(in_ch, out_ch, 3),
        nn.LeakyReLU(slope),
        nn.Conv2d(out_ch, out_ch, 3),
        nn.LeakyReLU(slope),
    )


def _center_crop(x, h, w):
    dy = (x.shape[-2] - h) // 2
    dx = (x.shape[-1] - w) // 2
    return x[..., dy:dy + h, dx:dx + w]


class UNet2D(nn.Module):

    def __init__(self, in_channels=3, n_classes=3, base_features=16, n_levels=7,
                 negative_slope=0.1, input_offset=-0.5):
        super().__init__()
        self.n_levels = n_levels
        self.input_offset = input_offset
        feats = [base_features * 2 ** i for i in range(n_levels + 1)]

        self.down = nn.ModuleList()
        ch = in_channels
        for f in feats:
            self.down.append(_conv_block(ch, f, negative_slope))
            ch = f

        self.upconv = nn.ModuleList()
        self.up = nn.ModuleList()
        for i in reversed(range(n_levels)):
            self.upconv.append(nn.ConvTranspose2d(feats[i + 1], feats[i], 2, stride=2))
            self.up.append(_conv_block(2 * feats[i], feats[i], negative_slope))
        self.act = nn.LeakyReLU(negative_slope)

        self.score = nn.Conv2d(feats[0], n_classes, 1)

    def forward(self, x):
        x = x + self.input_offset
        skips = []
        for i, block in enumerate(self.down):
            x = block(x)
            if i < self.n_levels:
                skips.append(x)
                x = F.max_pool2d(x, 2, 2)
        for upconv, block in zip(self.upconv, self.up):
            x = self.act(upconv(x))
            skip = skips.pop()
            x = block(torch.cat([x, _center_crop(skip, x.shape[-2], x.shape[-1])], 1))
        return self.score(x)

    # ------------------------------------------------------------------
    # Shape helpers (valid convolutions => output smaller than input)

    @property
    def downsample_factor(self):
        return 2 ** self.n_levels

    def output_shape(self, n):
        """Output size for a given input size, or None if n is not valid."""
        for _ in range(self.n_levels):
            n -= 4
            if n <= 0 or n % 2:
                return None
            n //= 2
        n -= 4
        if n <= 0:
            return None
        for _ in range(self.n_levels):
            n = 2 * n - 4
        return n if n > 0 else None

    def input_shape(self, out):
        """Input size needed to produce an output of size out."""
        n = out
        for _ in range(self.n_levels):
            n = (n + 4) // 2
        n += 4
        for _ in range(self.n_levels):
            n = 2 * n + 4
        return n

    def valid_output_shape(self, requested):
        """Largest valid output tile size <= requested (at least the minimum)."""
        smallest = next(o for o in map(self.output_shape, itertools.count(1)) if o)
        # Valid output sizes are smallest + k * downsample_factor
        k = max(0, (requested - smallest) // self.downsample_factor)
        return smallest + k * self.downsample_factor
