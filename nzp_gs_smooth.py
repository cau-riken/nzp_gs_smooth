# -*- coding: utf-8 -*-
"""
Apply Gauss-Seidel Iteration Scheme to smooth high-frequency distortions across a slice stack.

Nipype ready code - inner functions used inside the 'run' function.

Author: Alexander Woodward, Connectome Analysis Unit, RIKEN CBS, Wako, Japan
Email: alexander.woodward@riken.jp

"""

def run(in_dir, cur_dir, out_dir_images, out_dir_transforms,
        iterations, ants_thread_count):
    """Apply Gauss-Seidel Iteration Scheme algorithm to a folder of images in sequence.

    Args:
        in_dir: The folder of images to process. Images should be padded with leading zeroes and a starting index of 1, e.g. slice_001.tif
        cur_dir: The current directory to work from.
        out_dir_images: The directory to place the transformed images.
        out_dir_transforms: Directory for calculated transforms.
        iterations: Number of iterations of G.S. smoothing, 2,4 is recommended.
        ants_thread_count: Number of threads assigned to the antsRegistration program.

    """
    import numpy as np
    # TODO(A.W.): Remove dependency on OpenCV and use only SimpleITK calls
    import cv2
    import subprocess
    import SimpleITK as sitk
    import copy
    import os
    from nipype.interfaces.ants import Registration
    from natsort import natsorted

    def reg_run(fixed_image, moving_image, output_transform_prefix, output_warped_image, ants_thread_count):
        # os.environ['PATH']+=':/path_to_antsbin'
        reg = Registration()
        reg.inputs.fixed_image = fixed_image
        reg.inputs.moving_image = moving_image
        reg.inputs.output_transform_prefix = output_transform_prefix
        reg.inputs.transforms = ['SyN']
        reg.inputs.transform_parameters = [(0.01,)]
        reg.inputs.number_of_iterations = [[200, 200, 200, 200, 150, 50]]
        # reg.inputs.number_of_iterations = [[50,50,50,40,30,20]]
        reg.inputs.dimension = 2
        reg.inputs.num_threads = ants_thread_count
        reg.inputs.metric = ['Mattes']
        # Default (value ignored currently by ANTs)
        reg.inputs.metric_weight = [1]
        reg.inputs.radius_or_number_of_bins = [32]
        reg.inputs.sampling_strategy = ['Regular']
        reg.inputs.sampling_percentage = [1.0]  # 0.3]
        reg.inputs.convergence_threshold = [1.e-8]
        reg.inputs.convergence_window_size = [10]
        reg.inputs.smoothing_sigmas = [[6, 5, 4, 3, 2, 1]]
        reg.inputs.sigma_units = ['vox']
        reg.inputs.shrink_factors = [[6, 5, 4, 3, 2, 1]]
        reg.inputs.use_estimate_learning_rate_once = [True]
        reg.inputs.use_histogram_matching = [True]  # This is the default
        reg.inputs.output_warped_image = output_warped_image
        reg1 = copy.deepcopy(reg)
        reg1.cmdline
        reg1.run()

    # Be careful to scale by the pixel size when using both SimpleITK and OpenCV functions
    def deform_image(img_in, deform_in, scale_factor, out_filename):
        spacing = img_in.GetSpacing()
        a = sitk.GetArrayFromImage(img_in)
        a = a.astype(np.float32)
        c = sitk.GetArrayFromImage(deform_in)
        mx = np.zeros([deform_in.GetHeight(), deform_in.GetWidth()])
        my = np.zeros([deform_in.GetHeight(), deform_in.GetWidth()])
        my = my.astype(np.float32)
        mx = mx.astype(np.float32)
        for j in range(0, deform_in.GetHeight()):
            for i in range(0, deform_in.GetWidth()):
                pix = c[j, i]
                mx[j, i] = float(j)+scale_factor * \
                    float(pix[1]) * (1.0/spacing[1])
                my[j, i] = float(i)+scale_factor * \
                    float(pix[0]) * (1.0/spacing[0])

        out = cv2.remap(a, my, mx, cv2.INTER_CUBIC)
        out = sitk.GetImageFromArray(out)
        out.CopyInformation(img_in)
        sitk.WriteImage(out, out_filename)
        return out

    # copy from the input directory first
    names_in = subprocess.check_output(["ls "+in_dir], shell=True, text=True)
    names_in = names_in.split()
    # sort the strings using natural sort
    names_in = natsorted(names_in)
    img_count = len(names_in)
    
    # copy and duplicate first and last
    for i in range(0, img_count):
        #print("Attempting to load "+in_dir+names_in[i])
        img = sitk.ReadImage(in_dir+names_in[i])
        index = i+1
        sitk.WriteImage(
            img, cur_dir+"/input_with_boundary/slice_"+format(index, '04d')+".nii")
        if (index == 1):
            sitk.WriteImage(
                img, cur_dir+"/input_with_boundary/slice_0000.nii")
        elif (index == img_count):
            sitk.WriteImage(
                img, cur_dir+"/input_with_boundary/slice_"+format(index+1, '04d')+".nii")
    # Copy files into directories
    subprocess.call("cp "+cur_dir+"/input_with_boundary/* "+cur_dir+"/current_iter/", shell=True)
    print("Copied files to ./current_iter")
    subprocess.call("cp "+cur_dir+"/input_with_boundary/* "+cur_dir+"/prev_iter/", shell=True)
    print("Coped files to ./prev_iter")
    names_in = subprocess.check_output(["ls "+cur_dir+"/current_iter/"], shell=True, text=True)
    names_in = names_in.split()
    img_count = len(names_in)
    image_names = []

    # Neumann boundary already included
    for i in range(0, img_count):
        image_names.append(names_in[i])

    # Do nothing if iterations == 0
    if (iterations > 0):
        # iterations = 4
        Z = img_count
        for t in range(0, iterations):
            print("Starting iteration "+str(t+1)+" of "+str(iterations))
            if (t % 2 == 0):
                # increment
                for j in range(1, Z-1): # Z-1
                    # calculate transform between j-1 and j+1
                    # fixed is j-1 moving is j+1
                    reg_run(cur_dir+'/current_iter/'+image_names[j-1], cur_dir+'/current_iter/'+image_names[j+1],
                            cur_dir+'/registration_output_transform/u_', cur_dir+'/registration_output_image/output_warped_image.nii', ants_thread_count)
                    img_u = sitk.ReadImage(cur_dir+'/registration_output_transform/u_0Warp.nii.gz', sitk.sitkVectorFloat64)
                    # multiply it by half the deformation field
                    img_JM1 = sitk.ReadImage(cur_dir+'/current_iter/'+image_names[j+1])
                    # make sure to account for pixel scale factor
                    deform_image(img_JM1, img_u, 0.5, cur_dir+'/i_hat_image/output.nii')
                    # register Ij to IHat
                    reg_run(cur_dir+'/i_hat_image/output.nii', cur_dir+'/prev_iter/'+image_names[j],cur_dir+'/registration_output_transform/u_', cur_dir+'/registration_output_image/output_warped_image.nii', ants_thread_count)
                    #subprocess.check_output(['rm '+cur_dir+'/registration_output_transform/u_0Warp.nii'], shell=True)
                    # merge it with previous
                    img_uAccNew = sitk.ReadImage(cur_dir+'/registration_output_transform/u_0Warp.nii.gz', sitk.sitkVectorFloat64)
                    if (t == 0):
                        sitk.WriteImage(img_uAccNew, cur_dir+'/current_transforms/u'+str(j)+'_Warp.nii')
                    else:
                        img_uAcc = sitk.ReadImage(cur_dir+'/current_transforms/u'+str(j)+'_Warp.nii', sitk.sitkVectorFloat64)
                        img_uAccNew = img_uAcc+img_uAccNew
                        sitk.WriteImage(img_uAccNew, cur_dir+'/current_transforms/u'+str(j)+'_Warp.nii')
                    # update Ij using Ij0
                    img_JOrig = sitk.ReadImage(cur_dir+'/input_with_boundary/'+image_names[j])
                    deform_image(img_JOrig, img_uAccNew, 1.0, cur_dir+'/current_iter/'+image_names[j])
                # update boundaries
                img_b = sitk.ReadImage(cur_dir+'/current_iter/'+image_names[1])
                sitk.WriteImage(img_b, cur_dir+'/current_iter/'+image_names[0])
                img_b = sitk.ReadImage(cur_dir+'/current_iter/'+image_names[Z-2])
                sitk.WriteImage(img_b, cur_dir+'/current_iter/'+image_names[Z-1])
                # copy Input to InputPreviousIteration
                for j in range(0, Z):
                    img_b = sitk.ReadImage(cur_dir+'/current_iter/'+image_names[j])
                    sitk.WriteImage(img_b, cur_dir+'/prev_iter/'+image_names[j])
            else:
                # increment
                for j in range(Z-2, 0, -1):  # Z-1
                    # calculate transform between j-1 and j+1
                    # fixed is j-1 moving is j+1
                    reg_run(cur_dir+'/current_iter/'+image_names[j+1], cur_dir+'/current_iter/'+image_names[j-1],
                            cur_dir+'/registration_output_transform/u_', cur_dir+'/registration_output_image/output_warped_image.nii', ants_thread_count)
                    img_u = sitk.ReadImage(cur_dir+'/registration_output_transform/u_0Warp.nii', sitk.sitkVectorFloat64)
                    # multiply it by half the deformation field
                    img_JM1 = sitk.ReadImage(cur_dir+'/current_iter/'+image_names[j-1])
                    deform_image(img_JM1, img_u, 0.5, cur_dir+'/i_hat_image/output.nii')
                    # register Ij to IHat
                    reg_run(cur_dir+'/i_hat_image/output.nii', cur_dir+'/prev_iter/'+image_names[j],
                            cur_dir+'/registration_output_transform/u_', cur_dir+'/registration_output_image/output_warped_image.nii', ants_thread_count)
                    # merge it with previous
                    img_uAccNew = sitk.ReadImage(cur_dir+'/registration_output_transform/u_0Warp.nii.gz', sitk.sitkVectorFloat64)
                    if (t == 0):
                        sitk.WriteImage(img_uAccNew, cur_dir+'/current_transforms/u'+str(j)+'_Warp.nii')
                    else:
                        img_uAcc = sitk.ReadImage(cur_dir+'/current_transforms/u'+str(j)+'_Warp.nii', sitk.sitkVectorFloat64)
                        img_uAccNew = img_uAcc+img_uAccNew
                        sitk.WriteImage(img_uAccNew, cur_dir+'/current_transforms/u'+str(j)+'_Warp.nii')
                    # update Ij
                    img_JOrig = sitk.ReadImage(cur_dir+'/input_with_boundary/'+image_names[j])
                    deform_image(img_JOrig, img_uAccNew,1.0,cur_dir+'/current_iter/'+image_names[j])
                # update boundaries
                img_b = sitk.ReadImage(cur_dir+'/current_iter/'+image_names[1])
                sitk.WriteImage(img_b, cur_dir+'/current_iter/'+image_names[0])
                img_b = sitk.ReadImage(cur_dir+'/current_iter/'+image_names[Z-2])
                sitk.WriteImage(img_b, cur_dir+'/current_iter/'+image_names[Z-1])
                # copy Input to InputPreviousIteration
                for j in range(0, Z):
                    img_b = sitk.ReadImage(cur_dir+'/current_iter/'+image_names[j])
                    sitk.WriteImage(img_b, cur_dir+'/prev_iter/'+image_names[j])

    # copy results to output directories
    subprocess.call("cp "+cur_dir+"/current_iter/* "+out_dir_images, shell=True)
    # remove the first and last since they were boundary conditions
    subprocess.call("rm "+out_dir_images+"/slice_0000.nii", shell=True)
    subprocess.call("rm "+out_dir_images+"/slice_"+format(img_count-1, '04d')+".nii", shell=True)
    # copy the transforms to the correct output folder
    subprocess.call("cp "+cur_dir+"/current_transforms/* "+out_dir_transforms, shell=True)
    # Fix the indexing
    # copy from the input directory first
    names_in = subprocess.check_output(["ls "+out_dir_images], shell=True, text=True)
    names_in = names_in.split()
    img_count = len(names_in)
    # these images should be ordered based on the file naming convention
    for i in range(0, img_count):
        print("Attempting to load "+out_dir_images+"/"+names_in[i])
        img = sitk.ReadImage(out_dir_images+"/"+names_in[i])
        sitk.WriteImage(img, out_dir_images+"/out_"+format(i+1, '04d')+".nii")
    # finally remove the old images
    subprocess.call("rm "+out_dir_images+"/slice_*", shell=True)
    return out_dir_images
