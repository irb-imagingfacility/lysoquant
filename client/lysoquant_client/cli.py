"""lysoquant-segment: segment a TIFF with a LysoQuant server from the command line"""

import argparse
import sys

import numpy as np
import tifffile

from .client import LysoQuantError, check_server, make_rgb, segment


def _channel(stack, number):
    """Channel number (from 1) of a TIFF (c, y, x) or (z, c, y, x)"""
    return stack[..., number - 1, :, :]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image", nargs="?", help="TIFF with the channels as (c, y, x) (or an RGB TIFF with --rgb)")
    parser.add_argument("-o", "--output", help="label image to write (TIFF, uint16: 0 background, 1 empty, 2 loaded)")
    parser.add_argument("--lysosome", type=int, default=2, help="lysosome channel, from 1 (default 2)")
    parser.add_argument("--protein", type=int, default=3, help="protein channel, from 1 (default 3)")
    parser.add_argument("--rgb", action="store_true", help="the image is already RGB: protein red, lysosomes green")
    parser.add_argument("--pixel-size", type=float, required=False, help="pixel size in microns (needed to segment)")
    parser.add_argument("--host", help="SSH server (default: the server runs on this computer)")
    parser.add_argument("--port", type=int, default=22)
    parser.add_argument("--username")
    parser.add_argument("--key", help="private key file")
    parser.add_argument("--folder", help="process folder of the server (default: looked for in the usual places)")
    parser.add_argument("--test", action="store_true", help="only test the connection and show what the server reports")
    args = parser.parse_args(argv)

    server = dict(host=args.host, port=args.port, username=args.username, key_file=args.key, folder=args.folder,
                  local=args.host is None)
    try:
        if args.test:
            info = check_server(**server)
            print(f"server running on {info['device']}, tile size {info['tile_size']} px, folder {info['folder']}")
            return 0
        if not (args.image and args.output and args.pixel_size):
            parser.error("image, --output and --pixel-size are needed to segment")

        image = tifffile.imread(args.image)
        rgb = image if args.rgb else make_rgb(_channel(image, args.protein), _channel(image, args.lysosome))
        result = segment(rgb, args.pixel_size, progress=lambda done, total: print(
            f"\rsegmenting, tile {done}/{total}", end="", file=sys.stderr), **server)
        print(file=sys.stderr)
        tifffile.imwrite(args.output, result.labels.astype(np.uint16), imagej=True,
                         resolution=(1 / result.pixel_size_um[1], 1 / result.pixel_size_um[0]),
                         metadata={"unit": "um"})
        print(f"{args.output}: labels {result.labels.shape}, pixel size {result.pixel_size_um[0]:.4g} um")
        return 0
    except (LysoQuantError, ConnectionError, OSError, TimeoutError) as e:
        print(f"lysoquant-segment: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("canceled", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
