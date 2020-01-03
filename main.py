# -*- coding: utf-8 -*-
"""
Standalone Gauss-Seidel Iteration Scheme to smooth high-frequency distortions across a slice stack

Entry point.

Author: Alexander Woodward, Connectome Analysis Unit, RIKEN CBS, Wako, Japan
Email: alexander.woodward@riken.jp

"""

import os
import sys
import getopt
import json
import shutil
from datetime import datetime
import argparse
import nzp_gs_smooth

def main(inputFolder): 
    print("Creating required folders.")
    print("This removes existing folders and data. Be careful!")
    cwd = os.getcwd()
    working_dir = cwd + "/working/"
    out_images = cwd+'/out_images'
    out_transforms = cwd+'/out_transforms'
    if os.path.exists(working_dir) is True:
        shutil.rmtree(working_dir)
    if os.path.exists(out_images) is True:
        shutil.rmtree(out_images)
    if os.path.exists(out_transforms) is True:
        shutil.rmtree(out_transforms)
    os.mkdir(out_images)#final output images
    os.mkdir(out_transforms)#final output transforms
    os.mkdir(working_dir)
    os.mkdir(working_dir+'/input_with_boundary')#OriginalImagesWithB
    os.mkdir(working_dir+'/current_iter')#Input
    os.mkdir(working_dir+'/prev_iter')#InputPreviousIteration
    os.mkdir(working_dir+'/registration_output_transform')#OutputTransforms
    os.mkdir(working_dir+'/current_transforms')#OutputTransformsU
    os.mkdir(working_dir+'/i_hat_image')#OutputIHAT
    os.mkdir(working_dir+'/registration_output_image')#Output

    # multiple of 2 is best for 'iterations'
    nzp_gs_smooth.run(inputFolder, working_dir, out_images, out_transforms, 2, 20)

if __name__ == '__main__':
    now = datetime.now()
    current_time_start = now.strftime("%H:%M:%S")
    print("Current time = ", current_time_start)
    parser = argparse.ArgumentParser(description='Standalone version of gauss-seidel smoothing of image stack.')
    parser.add_argument('inputFolder', metavar='inputFolder',
                    help='a folder with prealigned images')
    args = parser.parse_args()
    main(args.inputFolder)
    now = datetime.now()
    current_time_end = now.strftime("%H:%M:%S")
    print("Start time = ", current_time_start)
    print("End time = ", current_time_end)
