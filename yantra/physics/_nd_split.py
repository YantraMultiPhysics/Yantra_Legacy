# -*- coding: utf-8 -*-
"""
@author: ggoussar
"""

import numpy as np
from copy import deepcopy as copy


def refs(iterable):
    if type(iterable) == dict:
        for key in iterable:
            yield key
    else:
        i = 0
        maxi = len(iterable)
        while i < maxi:
            yield i
            i += 1


def nd_blocksplit_pattern(dim, split):
    target_dim_per_part = list(map(lambda x, y: x/y, dim, split))
    numlarge = list(map(lambda x, y, z: x - y*z, dim, target_dim_per_part, split))
    size = [np.ones((x,)) for x in split]
    for d in range(len(dim)):
        size[d][:numlarge[d]] = target_dim_per_part[d]+1
        size[d][numlarge[d]:] = target_dim_per_part[d]
    return size


def nd_tuple(dims, nd = 3):
    result = [1]*nd
    try:
        for i in range(-1,-nd-1,-1):
            result[i] = dims[i]
    except IndexError:
        pass
    except TypeError:
        pass
    return tuple(result)


def block_dims(sizes, indices):
    nd = len(sizes)
    result = [1]*nd
    i = 0
    for index in indices:
        result[i] = sizes[i][index]
        i += 1
    return tuple(result)


def nd_indices(max_values):
    complete = False
    nd = len(max_values)
    current = [0]*nd
    while not complete:
        yield tuple(current)
        i = nd-1
        current[i] += 1
        while current[i] > max_values[i]-1:
            current[i] = 0
            i -= 1
            if i >= 0:
                current[i] += 1
            else:
                complete = True


def blocksplit_nd(data, split, dim = None):
    orig_dim = ()
    try:
        orig_dim = data.shape
    except AttributeError:
        return [data]*np.prod(split)
    
    if dim == None:
        dim = data.shape
    if orig_dim != dim:
        data = data.reshape(dim)
    nd = max(len(split),len(dim))
    dim = nd_tuple(dim, nd)
    split = nd_tuple(split, nd)
    sizes = nd_blocksplit_pattern(dim, split)
    result = dimsplit(data, sizes)
    if orig_dim != dim:
        dim_delta = len(dim) - len(orig_dim)
        if dim_delta != 0:
            for i in range(len(result)):
                fd = 1
                shape = result[i].shape
                for k in range(dim_delta):
                    fd *= shape[k]
                target_dim = (fd*shape[dim_delta],) + shape[dim_delta:]
                result[i] = result[i].reshape(target_dim)
        else:
            errmsg = "Cannot reasonably transform shape %s into shape %s"% (
                dim, orig_dim)
            raise RuntimeError(errmsg)   
    return result


def blockmerge_nd(target, parts, split, dim = None):
    orig_dim = ()
    try:
        orig_dim = target.shape
    except AttributeError:
        return target

    if dim == None:
        dim = target.shape
    if orig_dim != dim:
        target = target.reshape(dim)
    nd = max(len(split),len(dim))
    dim = nd_tuple(dim, nd)
    split = nd_tuple(split, nd)
    sizes = nd_blocksplit_pattern(dim, split)
    dimmerge(target, iter(parts), sizes)
    return target


def dimsplit(array, sizes, axis = 0):
    maxaxis = len(sizes)-1
    result = []
    actual_sizes = copy(sizes[axis])
    for i in range(1,len(actual_sizes)):
        actual_sizes[i] += actual_sizes[i-1]
    tmp = np.split(array, actual_sizes[:-1], axis)
    if axis < maxaxis:
        for part in tmp:
            result.extend(dimsplit(part, sizes, axis+1))
    else:
        result = tmp
    return result


def dimmerge(target, parts_it, sizes, axis = 0, sup_slice = ()):
    maxaxis = len(sizes)
    if axis < maxaxis:
        slice_begin = 0
        for index in range(len(sizes[axis])):
            slice_end = slice_begin + sizes[axis][index]
            sub_slice = sup_slice + (slice(slice_begin, slice_end), )
            dimmerge(target, parts_it, sizes, axis+1, sub_slice)
            slice_begin = slice_end
    else:
        part = next(parts_it)
        ss_shape = ()
        for i in range(len(sup_slice)):
            ss_shape += (sup_slice[i].stop - sup_slice[i].start, )
        if ss_shape != part.shape:
            part = part.reshape(ss_shape)
        target[sup_slice] = part


def nd_split(data, split, dim = None):
    parts = []
    if type(data) != np.ndarray:
        try:
            ncores = np.prod(split)
            parts = [None]*ncores
            if type(data) != dict:
                for proc in range(ncores):
                    parts[proc] = [None]*len(data) 
            else:
                for proc in range(ncores):
                    parts[proc] = {}
            for key in refs(data):
                values = data[key]
                part = nd_split(values, split, dim)
                for proc in range(ncores):
                    parts[proc][key] = part[proc]
        except TypeError:
            parts = blocksplit_nd(data, split, dim)
    else:
        parts = blocksplit_nd(data, split, dim)
    return parts


def nd_merge(data, parts, split, dim = None):
    result = data
    if type(data) != np.ndarray:
        try:
            for key in refs(data):
                target_parts = []
                for part in parts:
                    target_parts.append(part[key])
                result[key] = nd_merge(data[key], target_parts, split, dim)
        except TypeError:
            result = blockmerge_nd(data, parts, split, dim)
    else:
        result = blockmerge_nd(data, parts, split, dim)
    return result


# Use example
if __name__ == '__main__':
    a= np.array([[111,112],[121,122]])
    b= np.array([[211,212],[221,222]])
    c= np.array([[311,312],[321,322]])
    d= np.array([[411,412],[421,422]])

    vals = {'a':a,'b':b,'c':c,'d':d}
    abcd_parts = copy(nd_split(vals,(2,2),(2,2)))
    for i in abcd_parts:
        for key in refs(i):
            i[key] = i[key]%100+(i[key]/100)*1000
    vals = nd_merge(vals, abcd_parts,(2,2),(2,2))
    print(vals)
    #print nd_split([a,b,c,d], (2,2))
