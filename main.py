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

def create_dirs(): 
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
    os.mkdir(out_images)#Final output images
    os.mkdir(out_transforms)#Final output transforms
    os.mkdir(working_dir)
    os.mkdir(working_dir+'/input_with_boundary')
    os.mkdir(working_dir+'/current_iter')
    os.mkdir(working_dir+'/prev_iter')
    os.mkdir(working_dir+'/registration_output_transform')
    os.mkdir(working_dir+'/current_transforms')
    os.mkdir(working_dir+'/i_hat_image')
    os.mkdir(working_dir+'/registration_output_image')

    dirs = {
        'working_dir': working_dir,
        'out_images': out_images,
        'out_transforms': out_transforms
    }
    return dirs
    
if __name__ == '__main__':
    start_time = datetime.now()
    
    parser = argparse.ArgumentParser(description='Standalone version of gauss-seidel smoothing of image stack.')
    parser.add_argument('input_dir', metavar='input_dir',
                    help='A directory with prealigned images')
    args = parser.parse_args()
    dirs = create_dirs()

    # Run the procedure
    nzp_gs_smooth.run(args.input_dir, dirs['working_dir'], dirs['out_images'], dirs['out_transforms'], 2, 20)

    end_time = datetime.now()
    
    print("Total running time = ", end_time-start_time)
    
