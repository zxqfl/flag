#!/bin/bash
set -eux
cd "$(dirname $0)"

cargo run --release
open test.png

# ffmpeg -y -framerate 24 -pattern_type glob -i 'test-*.png' -vf "scale=iw/8:ih/8,reverse"  -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
# open out.mp4
