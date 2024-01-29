from PIL import Image, ImageTransform
from common import ImageLoader
import numpy as np
import glob
import os
import sys
import shutil

import matplotlib.pyplot as plt

import skimage.restoration

def downscale_mean(image, factor):
    image = image[:image.shape[0]//factor*factor,:image.shape[1]//factor*factor]
    image = image.reshape(image.shape[0]//factor, factor, image.shape[1]//factor, factor, *image.shape[2:])
    return np.mean(image, axis=(1,3)).astype(image.dtype)

class CompositeImage:
    def __init__(self):
        self.image = None
        self.x = None
        self.y = None
        self.width = None
        self.height = None
    
    def xres(self):
        if self.image is None:
            return 1
        return self.image.shape[0]/self.width

    def yres(self):
        if self.image is None:
            return 1
        return self.image.shape[1]/self.height

    def add_image(self, newimage, x, y, width=None, height=None):
        #print ('adding image', newimage.shape, x, y, width, height, self.x, self.y, self.width, self.height, self.xres(), self.yres())
        if width is None:
            width = newimage.shape[0]/self.xres()
        if height is None:
            height = newimage.shape[1]/self.yres()
        
        if (self.image is None):
            self.image = newimage
            self.mask = np.ones(newimage.shape, dtype='uint8')
            self.x, self.y, self.width, self.height = x, y, width, height
        else:
            if abs(newimage.shape[0]/width - self.xres()) > 1/self.xres() or abs(newimage.shape[1]/height - self.yres()) > 1/self.yres():
                print (newimage.shape[0]/width, self.xres(), newimage.shape[1]/height, self.yres())
                newimage, width, height = self.rescale(newimage, width, height)
            pix_x, pix_y = int((x - self.x) * self.xres()), int((y - self.y) * self.yres())

            if pix_x < 0 or pix_x+newimage.shape[0] > self.image.shape[0] or pix_y < 0 or pix_y+newimage.shape[1] > self.image.shape[1]:
                newshape = (max(self.image.shape[0] + max(0,-pix_x), pix_x + newimage.shape[0]), max(self.image.shape[1] + max(0,-pix_y), pix_y + newimage.shape[1]), *self.image.shape[2:])
                newrootimage = np.zeros(newshape, dtype=self.image.dtype)
                newrootmask = np.zeros(newshape, dtype='uint8')
                negpix_x, negpix_y = max(0,-pix_x), max(0,-pix_y)
                newrootimage[negpix_x:negpix_x+self.image.shape[0], negpix_y:negpix_y+self.image.shape[1]] = self.image
                newrootmask[negpix_x:negpix_x+self.image.shape[0], negpix_y:negpix_y+self.image.shape[1]] = self.mask
                
                self.x = min(self.x, x)
                self.y = min(self.y, y)
                self.width = newrootimage.shape[0]/self.xres()
                self.height = newrootimage.shape[1]/self.yres()
                
                self.image = newrootimage
                self.mask = newrootmask

                pix_x, pix_y = int((x - self.x) * self.xres()), int((y - self.y) * self.yres())
            
            x0,x1,y0,y1 = pix_x,pix_x+newimage.shape[0], pix_y,pix_y+newimage.shape[1]
            #mask = self.mask[x0:x1,y0:y1]
            #self.image[x0:x1,y0:y1][mask!=0] //= 2
            #self.image[x0:x1,y0:y1][mask!=0] = (self.image[x0:x1,y0:y1][mask!=0].astype(int) * mask[mask!=0] // (mask[mask!=0]+1)).astype('uint8')
            self.image[x0:x1,y0:y1] = newimage
            self.mask[x0:x1,y0:y1] = 1
            #self.image[pix_x:pix_x+newimage.shape[0], pix_y:pix_y+newimage.shape[1]] = newimage
            #self.image[pix_x:pix_x+newimage.shape[0], pix_y:pix_y+newimage.shape[1]] += newimage
            #self.mask[pix_x:pix_x+newimage.shape[0], pix_y:pix_y+newimage.shape[1]] += 1

    def full_image(self, fill_val=0):
        full_image = self.image.copy()
        full_image[self.mask==0] = fill_val
        return full_image
    
    def best_position(self, image, start_x, start_y, width=None, height=None, search_radius=25, grow_image=False):
        if grow_image:
            newimage = np.zeros((image.shape[0]+2, image.shape[1]+2, *image.shape[2:]), dtype=image.dtype)
            newimage[1:-1,1:-1] = image
            newimage[0,1:-1] = image[0,:]
            newimage[-1,1:-1] = image[-1,:]
            newimage[1:-1,0] = image[:,0]
            newimage[1:-1,-1] = image[:,-1]
            start_x -= 1/self.xres()
            start_y -= 1/self.yres()
            image = newimage
        
        if width is None:
            width = image.shape[0]/self.xres()
        if height is None:
            height = image.shape[1]/self.yres()
        
        start_pix_x, start_pix_y = int((start_x - self.x) * self.xres()), int((start_y - self.y) * self.yres())

        if type(search_radius) == int:
            search_radius = search_radius, search_radius
        
        best_score = 0
        best_pos = 0,0

        for xoff in range(-search_radius[0], search_radius[0]+1):
            for yoff in range(-search_radius[1], search_radius[1]+1):
                pix_x, pix_y = start_pix_x + xoff, start_pix_y + yoff
                
                img_x, img_y = max(0,-pix_x), max(0,-pix_y)
                pix_x, pix_y = max(0,pix_x), max(0,pix_y)
                width, height = max(0,image.shape[0]-img_x), max(0,image.shape[1]-img_y)
                #width = min(self.image.shape[0], pix_x+width) - self.image.shape[0] + width
                #height = min(self.image.shape[1], pix_y+height) - self.image.shape[1] + height
                
                maskslice = self.mask[pix_x:pix_x+width, pix_y:pix_y+height]
                rootslice = self.image[pix_x:pix_x+width, pix_y:pix_y+height]
                width, height = rootslice.shape[:2]
                imageslice = image[img_x:img_x+width, img_y:img_y+height]
                #print (imageslice.shape, maskslice.shape, rootslice.shape, image.shape, width, height, img_x, img_y, pix_x, pix_y)
                if np.sum(maskslice) > 0:
                    score = np.sum(rootslice[maskslice] * imageslice[maskslice]) / np.sum(maskslice)
                    if score > best_score:
                        best_score = score
                        best_pos = xoff, yoff
        
        if grow_image:
            return (start_pix_x + best_pos[0] + 1) / self.xres(), (start_pix_y + best_pos[1] + 1) / self.yres()
        return best_pos[0] / self.xres(), best_pos[1] / self.yres()





def illum_correction(images):
    images = np.array(images)
    background = np.min(images, axis=0)
    #signal = np.percentile(images, axis=0).astype(images[0].dtype) - background
    images = (images - background)
    signal_val = np.percentile(images, 99.9, axis=(0,1,2))
    images = images / signal_val
    images[images<0] = 0
    images[images>1] = 1
    return images

def stitch_grid_manual(images, param_list, grid_size, x_skew, y_skew):
    composite = CompositeImage()

    for image, params in zip(images, param_list):
        image = image[2:-2,...]
        i = params['section']-1
        x = i//grid_size
        y = i%grid_size
        x, y = x*image.shape[0] + y*x_skew, y*image.shape[1] + x*y_skew
        composite.add_image(image, x, y)

    return composite


def stitch_grid_metadata(images, param_list, imageloader, metadata_func):
    composite = CompositeImage()

    for image, params in zip(images, param_list):
        position = [*metadata_func(image, imageloader.to_path(params))]
        #if (not composite.image is None):
            #best_position = composite.best_position(image, *position, search_radius=20)
            #print (best_position, "best pos")
            #position[0] += best_position[0]
            #position[1] += best_position[1]
        composite.add_image(image, *position)

    return composite

def stitch_grid_m2stitch(images, param_list, imageloader, grid_pos_func, metadata_func = None, silent=False, **kwargs):
    if silent:
        import tqdm
        real_tqdm = tqdm.tqdm
        def replacement(iterable, *args, **kwargs):
            return iterable
        tqdm.tqdm = replacement
    import m2stitch
    composite = CompositeImage()
    
    stitch_images = np.array(images)
    if len(stitch_images.shape) > 3:
        stitch_images = np.sum(stitch_images, axis=tuple(range(3,len(stitch_images.shape))))
    
    positions = np.array([grid_pos_func(images[i], param_list[i], imageloader.to_path(param_list[i]))[:2] for i in range(len(images))])
    positions -= positions.min(axis=0)

    #if os.path.exists('result.csv'):
    #    result = pandas.from_csv('result.csv')
    #else:
    if True:
        if metadata_func:
            guesses = np.array([metadata_func(images[i], imageloader.to_path(param_list[i])) for i in range(len(images))])
            result, result_dict = m2stitch.stitch_images(stitch_images, position_indices=positions, position_initial_guesses=guesses, **kwargs)
        else:
            result, result_dict = m2stitch.stitch_images(stitch_images, position_indices=positions, **kwargs)
        #result.to_csv('result.csv')

    for i in range(len(images)):
        composite.add_image(images[i], result.x_pos[i], result.y_pos[i], images[i].shape[0], images[i].shape[1])

    if silent:
        tqdm.tqdm = real_tqdm
    
    return composite

def stitch_test():
    imageloader = ImageLoader()
    imageloader.image_path = 'tmp/'
    imageloader.file_pattern = 'slide{x:d}x{y:d}.png'
    imageloader.load()
    
    images, params = imageloader.load_params()
    
    rows = [param['x'] for param in params]
    cols = [param['y'] for param in params]
    
    positions = [(param['y'], param['x']) for param in params]
    guesses = [(param['x']*100, param['y']*100) for param in params]

    for i in range(len(rows)):
        width, height = images[i].shape[:2]
        #image = np.stack([images[i]]*8, axis=1)
        #image = np.stack([image]*8, axis=3)
        #images[i] = image.reshape(width*8, height*8, -1)
        #guesses.append((rows[i]*100, cols[i]*100))
        #fig, axis = plt.subplots()
        #axis.imshow(images[i].sum(axis=2))
        #axis.set_title("{}x{}".format(rows[i], cols[i]))
    #plt.show()
    
    import m2stitch
    stitch_images = np.array(images).sum(axis=3).astype('uint16')
    print (stitch_images.max(), stitch_images.min(), stitch_images.dtype)
    print (stitch_images.shape)
    print (positions)
    print (guesses)
    result, result_dict = m2stitch.stitch_images(stitch_images, position_indices=positions, ncc_threshold=0.1)#, position_initial_guess=guesses)

    print (result)

    composite = CompositeImage()
    for i in range(len(images)):
        composite.add_image(images[i], result.x_pos[i], result.y_pos[i], images[i].shape[0], images[i].shape[1])
    
    skimage.io.imsave('test_stitch.png', composite.full_image())


