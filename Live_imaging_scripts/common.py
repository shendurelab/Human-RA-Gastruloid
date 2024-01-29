import numpy as np
import matplotlib.pyplot as plt
#from cellpose import models
#from cellpose.io import imread
import sys
import glob
import re
import parse
import os
import time
import shutil
import skimage.io
import filecmp
import pickle
import atexit

#from stardist.models import StarDist2D
#from csbdeep.utils import normalize

#image_path = "images/"
#image_path = "/mnt/nas2/Timecourses/6-8-22_A549_tripleKI_fluorbrite_drugs_v3/6-8-22_A549_tripleKI_fluorbrite_drugs_v3_"
#image_path = "/mnt/nas2/Timecourses/6-24-22_A549_tripleKIC1_drugs_v5/6-24-22_A549_tripleKIC1_drugs_v5_"

class ImageLoader:
    def __init__(self, pattern=None):#'t{time:d}_{strength}0{drug:d}_s{section:d}_w{layer:d}.TIF'):
        #self.image_path = 'images/'
        self.pattern = pattern
        if self.pattern is not None:
            self.load()
        #self.file_pattern = 't{time:d}_{strength}0{drug:d}_s{section:d}_w{layer:d}.TIF'

    def load(self):
        wildcard_pattern = re.sub('\{[^\}]*\}', '*', self.pattern)
        print (wildcard_pattern)
        self.paths = sorted(glob.glob(wildcard_pattern))
        #all_paths = [os.path.join(dp, f) for dp, dn, fn in os.walk(self.image_path) for f in fn]
        #self.paths = sorted([path for path in all_paths if self.to_params(path) != None])
        #class Default(dict):
        #    def __missing__(self, key):
        #        return '*'
        #self.paths = sorted(glob.glob(self.to_path(Default())))
        self.param_ranges = {}
        for path in self.paths:
            params = self.to_params(path)
            for name in params:
                if name not in self.param_ranges:
                    self.param_ranges[name] = [params[name]]
                elif params[name] not in self.param_ranges[name]:
                    self.param_ranges[name].append(params[name])

        for name in self.param_ranges:
            self.param_ranges[name].sort()
    
    def to_path(self, params=None, **kwargs):
        #class RemoveDict(dict):
        #    def __getitem__(self, key):
        #        return self.pop(key)
        #params = RemoveDict(params.copy())
        #path = file_pattern.format_map(params)
        #print (path, params)
        #path = "t{}_{}0{}_s{}_w{}".format(params.pop('time'), params.pop('strength'), params.pop('drug'), params.pop('section'), params.pop('layer'))
        #path += '_'.join([name+str(val) for name,val in params])
        #path = self.image_path + self.file_pattern.format_map(params)
        path = self.pattern.format_map(params or kwargs)
        return path

    def make_path(self, params=None, **kwargs):
        path = self.to_path(params or kwargs)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        return path
    
    def to_params(self, path):
        result = parse.parse(self.pattern, path, case_sensitive=True)
        #result = parse.parse(self.file_pattern, path.replace(self.image_path, ''), case_sensitive=True)
        #param_words = path.replace(self.image_path,'').split('.')[0].split('_')
        #params = dict(zip("time drug section layer".split(), [int(str[1:]) for str in param_words]))
        #params['strength'] = param_words[1][0]
        return None if result == None else result.named
    
    def all_params(self, **params):
        all_params = []
        for name in self.param_ranges:
            if name not in params:
                for val in self.param_ranges[name]:
                    new_params = params.copy()
                    new_params[name] = val
                    all_params += self.all_params(**new_params)
                break
        else:
            for name in params:
                if type(params[name]) != str and (hasattr(params[name], '__iter__') or hasattr(params[name], '__getitem__')):
                    values = params[name]
                    for val in values:
                        params[name] = val
                        all_params += self.all_params(**params)
                    break
            else:
                return [params]
        return all_params

    def load_params(self, **params):
        params = self.all_params(**params)
        images = [skimage.io.imread(self.to_path(param)) for param in params]
        return images, params

    def load_param_stack(self, stack_param=None, **params):
        params = self.all_params(**params)
        images = [skimage.io.imread(self.to_path(param)) for param in params]
        if stack_param != None:
            if stack_param in self.param_ranges and stack_param not in params:
                new_images_list = []
                for cur_param in self.param_ranges[stack_param]:
                    new_images_list.append([image for image,param in zip(images,params) if param[stack_param] == cur_param])
                images = [np.stack(image_list, axis=-1) for image_list in zip(*new_images_list)]
                params = [param for param in params if param[stack_param] == self.param_ranges[stack_param][0]]
        return images, params


def backup_folder(basedir):
    history_dir = basedir + '/history/'
    timestamp = time.strftime('%Y-%m-%d-{:04}')
    for i in range(10000):
        if not os.path.exists(history_dir + timestamp.format(i)):
            timestamp = timestamp.format(i)
            curdir = history_dir + timestamp
            break

    all_files = [filename for filename in os.listdir(basedir) if filename != 'history']

    if (os.path.exists(history_dir) and len(os.listdir(history_dir))):
        lastdir = [dir for dir in sorted(os.listdir(history_dir)) if dir not in ['byname', 'latest']]
        lastdir = history_dir + lastdir[-1]
        for filename in all_files:
            lastfile, basefile, curfile = [dir + '/' + filename for dir in [lastdir, basedir, curdir]]
            if not os.path.exists(lastfile) or not filecmp.cmp(basefile, lastfile):
            #if not os.path.exists(lastfile) or os.path.getmtime(basefile) > os.path.getmtime(lastfile):
                os.makedirs(os.path.dirname(curfile), exist_ok=True)
                shutil.copyfile(basefile, curfile)
                byname = history_dir + 'byname/' + filename + '/' + timestamp + '-' + filename
                os.makedirs(os.path.dirname(byname), exist_ok=True)
                os.link(curfile, byname)
    elif len(all_files):
        os.makedirs(curdir, exist_ok=True)
        for filename in all_files:
            shutil.copyfile(basedir + '/' + filename, curdir + '/' + filename)
            byname = history_dir + 'byname/' + filename + '/' + timestamp + '-' + filename
            os.makedirs(os.path.dirname(byname), exist_ok=True)
            os.link(curdir + '/' + filename, byname)

def mark_backup(basedir):
    atexit.register(backup_folder, basedir)

def load(name, calc_func=lambda: None):
    if os.path.exists(name):
        with open(name, 'rb') as ifile:
            return pickle.load(ifile)
    obj = calc_func()
    with open(name, 'wb') as ofile:
        pickle.dump(obj, ofile)
    return obj

