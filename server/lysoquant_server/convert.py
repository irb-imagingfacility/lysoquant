"""
Convert a caffe_unet model (.modeldef.h5 + .caffemodel.h5) into a single
PyTorch file that can be used by lysoquant-predict.

Usage:
    lysoquant-convert 2D-7-16.modeldef.h5 2D-7-16.caffemodel.h5 lysoquant.pt
"""

import argparse
import re

import numpy as np
import torch

from .unet import UNet2D

FORMAT = "lysoquant-unet-v1"


def _str(value):
    if isinstance(value, bytes):
        return value.decode("utf-8").rstrip("\x00")
    return str(value).rstrip("\x00")


def read_modeldef(path):
    import h5py      # only needed for the conversion, not for serving

    with h5py.File(path, "r") as f:
        p = f["unet_param"]
        prototxt = _str(f["model_prototxt"][()])
        config = {
            "name": _str(p["name"][()]),
            "description": _str(p["description"][()]),
            "class_names": [_str(c) for c in p["classNames"][()]],
            "element_size_um": [float(v) for v in p["element_size_um"][()]],
            "normalization_type": int(p["normalization_type"][()]),
            "padding": _str(p["padding"][()]),
        }

    slopes = set(float(s) for s in re.findall(r"negative_slope:\s*([-\d.eE]+)", prototxt))
    if len(slopes) != 1:
        raise ValueError(f"Expected a single LeakyReLU slope, found {slopes}")
    offset = re.search(r'type:\s*"ValueTransformation".*?offset\s*{\s*v:\s*([-\d.eE]+)', prototxt, re.S)
    config["negative_slope"] = slopes.pop()
    config["input_offset"] = float(offset.group(1)) if offset else 0.0
    return config


def read_caffe_weights(path, n_levels):
    """Map caffe_unet layer names to UNet2D parameter names."""
    import h5py

    mapping = {}
    for i in range(n_levels + 1):
        mapping[f"conv_d{i}a-b"] = f"down.{i}.0"
        mapping[f"conv_d{i}b-c"] = f"down.{i}.2"
    for i in range(n_levels):
        j = n_levels - 1 - i  # decoder modules are stored from the bottom up
        src = f"d{n_levels}c" if i == n_levels - 1 else f"u{i + 1}d"
        mapping[f"upconv_{src}_u{i}a"] = f"upconv.{j}"
        mapping[f"conv_u{i}b-c"] = f"up.{j}.0"
        mapping[f"conv_u{i}c-d"] = f"up.{j}.2"
    mapping["conv_u0d-score"] = "score"

    state = {}
    with h5py.File(path, "r") as f:
        data = f["data"]
        for caffe_name, torch_name in mapping.items():
            # Caffe stores weights as /data/<layer>/0 and biases as /data/<layer>/1.
            # Both Convolution and Deconvolution layouts match PyTorch directly.
            state[f"{torch_name}.weight"] = torch.from_numpy(np.asarray(data[caffe_name]["0"], dtype=np.float32))
            state[f"{torch_name}.bias"] = torch.from_numpy(np.asarray(data[caffe_name]["1"], dtype=np.float32))
    return state


def count_levels(path):
    import h5py

    with h5py.File(path, "r") as f:
        names = list(f["data"].keys())
    return max(int(m.group(1)) for m in (re.match(r"conv_d(\d+)a-b$", n) for n in names) if m)


def convert(modeldef, caffemodel, output):
    config = read_modeldef(modeldef)
    n_levels = count_levels(caffemodel)
    state = read_caffe_weights(caffemodel, n_levels)

    config.update({
        "n_levels": n_levels,
        "in_channels": state["down.0.0.weight"].shape[1],
        "base_features": state["down.0.0.weight"].shape[0],
        "n_classes": state["score.weight"].shape[0],
    })

    net = build_model(config)
    net.load_state_dict(state, strict=True)  # fails loudly on any mismatch

    torch.save({"format": FORMAT, "config": config, "state_dict": net.state_dict()}, output)
    return config


def build_model(config):
    return UNet2D(
        in_channels=config["in_channels"],
        n_classes=config["n_classes"],
        base_features=config["base_features"],
        n_levels=config["n_levels"],
        negative_slope=config["negative_slope"],
        input_offset=config["input_offset"],
    )


def load_model(path, device="cpu"):
    try:
        checkpoint = torch.load(path, map_location=device, weights_only=True)
    except Exception:
        checkpoint = {}
    if not isinstance(checkpoint, dict) or checkpoint.get("format") != FORMAT:
        raise ValueError(f"{path} is not a LysoQuant PyTorch model (use lysoquant-convert)")
    net = build_model(checkpoint["config"])
    net.load_state_dict(checkpoint["state_dict"])
    return net.to(device).eval(), checkpoint["config"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("modeldef", help="caffe_unet model definition (.modeldef.h5)")
    parser.add_argument("caffemodel", help="caffe_unet weights (.caffemodel.h5)")
    parser.add_argument("output", help="output PyTorch file (.pt)")
    args = parser.parse_args()

    config = convert(args.modeldef, args.caffemodel, args.output)
    print(f"Converted '{config['name']}' ({config['n_levels']} levels, "
          f"{config['base_features']} base features, classes {config['class_names']}) -> {args.output}")


if __name__ == "__main__":
    main()
