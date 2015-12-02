#!/bin/bash

# $1 is the folder where the PNGS are found

#ffmpeg -qp 0 -framerate 30 -i test%04d.png -c:v mpeg4 -r 30 $1_pos.mp4
# ffmpeg -framerate 30 -i $1/test_pos1%04d.png -q:v 0  -c:v mpeg4 -r 30 $1_pos1.mp4
# ffmpeg -framerate 30 -i $1/test_pos2%04d.png -q:v 0  -c:v mpeg4 -r 30 $1_pos2.mp4
ffmpeg -framerate 30 -i $1/fig%04d.png -q:v 0 -c:v mpeg4 -r 30 $1.mp4
# rm -f test*.png



# -b 5000k to set bitrate of file.
# -q:v 0 to set quality? This is the shit!!!
# -crf set quality   23 is default 18 is almost lossless 51 is worst