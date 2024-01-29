from stitch_images import CompositeImage, stitch_grid_m2stitch, stitch_grid_manual, downscale_mean, illum_correction
from common import ImageLoader
import sys
import skimage.io
import numpy as np
import concurrent.futures
import os
import matplotlib.pyplot as plt

def get_position_leica(image, path):
    fulldata = open(path, 'rb').read()
    text = fulldata.partition(b'<MetaData>')[2].partition(b'</MetaData>')[0].decode()
    xres = float(text.partition('<prop id="spatial-calibration-x" type="float" value="')[2].partition('"/>')[0])
    yres = float(text.partition('<prop id="spatial-calibration-y" type="float" value="')[2].partition('"/>')[0])
    #pix_width = int(text.partition('<prop id="pixel-size-x" type="int" value="')[2].partition('"/>')[0])
    #pix_height = int(text.partition('<prop id="pixel-size-y" type="int" value="')[2].partition('"/>')[0])
    x = float(text.partition('<custom-prop id="ASI X" type="float" value="')[2].partition('"/>')[0])
    y = float(text.partition('<custom-prop id="ASI Y" type="float" value="')[2].partition('"/>')[0])
    return y, x, image.shape[0] * yres, image.shape[1] * xres


def stitch_leica_all(dirpath, grid_size=None, stack_param='layer'):
    """ Stitches all of the images in an experiment directory which is formatted from the leica microscope
    The supplied dirpath should be the path of the base directory which includes all the files.
    """
    imageloader = ImageLoader(dirpath + "{prefix}_t{time:d}_{strength}{drug:02d}_s{section:d}_w{layer:d}.TIF")
    imageloader.load()

    print (imageloader.param_ranges)
    
    #timecourse_path = os.path.basename(dirpath)
    timecourse_path = imageloader.param_ranges['prefix'][0]
    date = timecourse_path.partition('_')[0]

    if not grid_size:
        grid_size = len(imageloader.param_ranges['section'])**0.5
        assert int(grid_size) == grid_size, "Unable to estimate grid size, specify it as an argument"
        grid_size = int(grid_size)
    
    imagesaver = ImageLoader("rgbimages/" + timecourse_path + "/{resolution}/{strength}{drug:02d}/t{time:03d}_w{layer:d}.tif")
    imagesaver.load()

    stitch_param = 'section'
    
    if not stack_param is None:
        for params in imageloader.all_params(**{stitch_param: imageloader.param_ranges[stitch_param][0], stack_param: imageloader.param_ranges[stack_param][0]}):
            params.pop(stack_param)
            params.pop(stitch_param)
            stitch_leica(imageloader, imagesaver, grid_size, params, stack_param=stack_param)
    else:
        for params in imageloader.all_params(**{stitch_param: imageloader.param_ranges[stitch_param][0]}):
            params.pop(stitch_param)
            stitch_leica(imageloader, imagesaver, grid_size, params)

def expand_image(image, radius):
    width, height = image.shape[:2]
    x0, x1 = radius, radius+width
    y0, y1 = radius, radius+height
    new_image = np.zeros((width+radius*2, height+radius*2) + image.shape[2:], dtype=image.dtype)
    new_image[x0:x1, y0:y1] = image
    new_image[:x0, :y0] = image[:1, :1]
    new_image[x1:, :y0] = image[-1:, :1]
    new_image[:x0, y1:] = image[:1, -1:]
    new_image[x1:, y1:] = image[-1:, -1:]
    new_image[:x0, y0:y1] = image[:1, :]
    new_image[x1:, y0:y1] = image[-1:, :]
    new_image[x0:x1, :y0] = image[:, :1]
    new_image[x0:x1, y1:] = image[:, -1:]
    return new_image

def stitch_leica(imageloader, imagesaver, grid_size, static_params, stack_param=None, do_illum_correction=False, manual_adjust=False):
    """ Stitches a field of view of a timeseries microscopy experiment.
    The alignment is done by hand, if you want to visualize the alignment and tweak it, run with manual_adjust=True.
    imageloader is the Imageloader that supplies the individual tiles and imagesaver is the one
    that specifies the path where the output is saved.
    """
    params = static_params.copy()
    params['resolution'] = 'full'
    if not stack_param is None: params[stack_param] = imageloader.param_ranges[stack_param][0]
    print (os.path.exists(imagesaver.to_path(params)), imagesaver.to_path(params))
    #if (os.path.exists(imagesaver.to_path(params))):
        #return
    print ("starting stitching ", static_params)
    try:
        if not stack_param is None:
            images, param_list = imageloader.load_params(stack_param=stack_param, **static_params)
        else:
            images, param_list = imageloader.load_params(**static_params)
    except FileNotFoundError:
        return
    print (np.min(images, axis=(0,1,2)), np.percentile(images, 1, axis=(0,1,2)), np.percentile(images, 99, axis=(0,1,2)), np.max(images, axis=(0,1,2)))
    images = np.array(images)
    print (np.min(images), np.max(images))
    print ("VALS")
    if do_illum_correction:
        images = (illum_correction(images) * 255).astype('uint8')
    print (images.shape)

    min_val, max_val = np.percentile(images, (0.1, 99.9), axis=(0,1,2))
    
    x_skew = -45
    y_skew = 45
    looping = manual_adjust
    composite = stitch_grid_manual(images, param_list, grid_size, x_skew=x_skew, y_skew=y_skew)
    image = composite.image

    while looping:
        expand_radius = int(images[0].shape[0]*0.3) 
        expanded_images = np.array([expand_image(image, expand_radius) for image in images])
        #composite = stitch_grid_manual(images, param_list, grid_size, x_skew=x_skew, y_skew=y_skew)
        #composite = stitch_grid_metadata(images, param_list, imageloader, get_position_leica)
        def grid_pos(image, params, path):
            i = params['section']-1
            x = i//grid_size
            y = i%grid_size
            return x, y, 1, 1
        #composite = stitch_grid_m2stitch(expanded_images, param_list, imageloader, grid_pos, ncc_threshold=0.01)
        import m2stitch
        composite = CompositeImage()

        positions = np.array([grid_pos(images[i], param_list[i], imageloader.to_path(param_list[i]))[:2] for i in range(len(images))])
        positions -= positions.min(axis=0)

        result, result_dict = m2stitch.stitch_images(expanded_images, position_indices=positions, ncc_threshold=0.01)

        new_poses = zip(result.x_pos, result.y_pos)
        for i in range(len(images)):
            x, y = result.x_pos[i], result.y_pos[i]
            #x -= expand_radius * positions[i,0]
            #y -= expand_radius * positions[i,1]
            composite.add_image(expanded_images[i], x, y, expanded_images[i].shape[0], expanded_images[i].shape[1])
            #composite.add_image(images[i], result.x_pos[i], result.y_pos[i], images[i].shape[0], images[i].shape[1])

        image = composite.image
        plt.imshow(image)
        plt.show()
        inp = input('change: ')
        looping = len(inp)
        if looping:
            xdif, ydif = eval(inp)
            x_skew += xdif
            y_skew += ydif
    #params = param_list[0].copy()
    #params['resolution'] = 'x3'
    #skimage.io.imsave(imagesaver.make_path(params), downscale_mean(image, 3))
    params['resolution'] = 'full'
    skimage.io.imsave(imagesaver.make_path(params), image, compression='deflate')
    #skimage.io.imsave('tmp.png', image.astype('uint8'))

"""
This program stitches together images from the leica microscope. It should work
with the standard file format the microscope exports to.
Arguments:
    path to directory of files
    optional, size of grid.
"""

if __name__ == "__main__":
    assert len(sys.argv) > 1, "Expected first argument to be a path to the directory with images"

    stitch_dir = sys.argv[1]
    grid_size = None if len(sys.argv) < 3 else int(sys.argv[2])
    stitch_leica_all(stitch_dir, grid_size=grid_size, stack_param=None)

