#!/usr/bin/python
# -*- coding: utf-8 -*-
#=======================================================================================
#This File is part of Yantra: A lattice Boltzmann method based tool for multiscale/
#multiphyics simulations
#=======================================================================================
#
#Copyright (C) 2016-2017  <Author> Ravi A. Patel <Email> ravee.a.patel@gmail.com
#
#This program is free software: you can redistribute it and/or modify it under the
#terms of the GNU General Public License as published by the Free Software 
#Foundation, either version 3 of the License, or any later version.
#This program is distributed in the hope that it will be useful, but WITHOUT ANY 
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
#PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#=======================================================================================

import numpy as np
import multiprocessing as mp

import time
__all__ =['list_split','list_merge','x_split','x_merge','y_split','y_merge','z_split',
          'z_merge','box_split','box_merge','run_parallel','get_merge_type']

          
def list_split(l,threads):
    """
    splits list into chunks of equal to number of threads
    
    Arguments
    ---------
    l: list
        list to be splitted
    threads: integer
        number of threads 
    
    Returns
    -------
    list
        list containing splitted chunks
    """
    n=int(np.ceil(len(l)/threads)) 
    out = [l[i:min(len(l),i +n)] for i in range(0, len(l), n)]
    if len(out)< threads:
        out.append([]*(threads-len(out)+1))
    return out


def list_merge(l):
    """
    merges list splitted by list_split
    
    Arguments
    ---------
    l: list
        list containing chunked list obtained from list_split
        
    Returns
    -------
    list
        merged list splitted using list_split
    """
    return [item for sublist in l for item in sublist]
    
def x_split(ndarray,threads):
    """
    splits arrays in number of threads specified along x direction
    
    Parameters
    ----------
    ndarray: ndarray
        array to be splitted
    threads: integer
        number of threads 
    
    Returns
    -------
    list of ndarray
        containing ndarray's chunks splitted into threads    
    """
    d = len(np.shape(ndarray))
    if d == 3:
        return np.array_split(ndarray,threads,2)
    elif d == 2:
        return np.array_split(ndarray,threads,1)
    else:
        raise ValueError("only 2D or 3D arrays can be splitted in x direction")
        
def y_split(ndarray,threads):
    """
    splits arrays in number of threads specified along y direction
    
    Parameters
    ----------
    ndarray: ndarray
        array to be splitted
    threads: integer
        number of threads 
    
    Returns
    -------
    list of ndarray
        containing ndarray's chunks splitted into threads    
    """
    d = len(np.shape(ndarray))
    if d == 3:
        return np.array_split(ndarray,threads,1)
    elif d == 2:
        return np.array_split(ndarray,threads,0)
    else:
        raise ValueError("only 2D or 3D arrays can be splitted in y direction")

def z_split(ndarray,threads):
    """
    splits arrays in number of threads specified along z direction
    
    Parameters
    ----------
    ndarray: ndarray
        array to be splitted
    threads: integer
        number of threads 
    
    Returns
    -------
    list of ndarray
        containing ndarray's chunks splitted into threads    
    """
    d = len(np.shape(ndarray))
    if d == 3:
        return np.array_split(ndarray,threads,0)
    else:
        raise ValueError("only 3D arrays can be splitted in z direction")

def box_split(ndarray,threads):
    """
    splits arrays in number of threads in cubic(square) chunks
    
    Parameters
    ----------
    ndarray: ndarray
        array to be splitted
    threads: integer
        number of threads 
    
    Returns
    -------
    list of ndarray
        containing ndarray's chunks splitted into threads    
    """
    d = len(np.shape(ndarray))
    if d == 2:
        chunksize = int(threads**0.5)
        return [y_split(i,chunksize) for i in x_split(ndarray,chunksize)]
    elif d == 3:
        pass
    
def x_merge(array_list):
    """
    merges arrays splitted in x direction using x_split
    
    Parameters
    ----------
    array_list: list of ndarray
        list containing ndarrays splitted using x_split
    
    Returns
    -------
    ndarray
        merged ndarray    
    """
    d = len(np.shape(array_list[0])) 
    if d == 3:
        return np.concatenate(array_list,2)
    elif d == 2:
        return np.concatenate(array_list,1)
    else:
        raise ValueError("only 2D or 3D arrays can be merged in x direction")

def y_merge(array_list):
    """
    merges arrays splitted in y direction using y_split
    
    Parameters
    ----------
    array_list: list of ndarray
        list containing ndarrays splitted using y_split
    
    Returns
    -------
    ndarray
        merged ndarray    
    """
    d = len(np.shape(array_list[0])) 
    if d == 3:
        return np.concatenate(array_list,1)
    elif d == 2:
        return np.concatenate(array_list,0)
    else:
        raise ValueError("only 2D or 3D arrays can be merged in y direction")

def z_merge(array_list):
    """
    merges arrays splitted in z direction using z_split
    
    Parameters
    ----------
    array_list: list of ndarray
        list containing ndarrays splitted using z_split
    
    Returns
    -------
    ndarray
        merged ndarray    
    """
    d = len(np.shape(array_list[0])) 
    if d == 3:
        return np.concatenate(array_list,0)
    else:
        raise ValueError("only 3D arrays can be merged in z direction")

def box_merge(array_list):
    """
    merges arrays splitted in cubic (square) chunks using box_split
    
    Parameters
    ----------
    array_list: list of ndarray
        list containing ndarrays splitted using box_split
    
    Returns
    -------
    ndarray
        merged ndarray    
    """
    d = len(np.shape(array_list[0][0])) 
    if d == 2:
        return x_merge([y_merge(i)for i in array_list])
    elif d == 3:
        pass


def get_merge_type(split_type):
    return split_type.replace('split','merge')

def parallel_func(func,pid,args):
    """
    returns output of func in list form with first element and pid (process id) and second element as output of the function
    
    Parameters:
    ----------
    func: function
        a function to be run
    pid: integer
         process id
    args: list
        arguments of the function. Arugment only functions are allowed at the moment.
    
    Returns:
    -------
        list 
    """
    out = func(*args)
    return [pid,out]
    
def run_parallel(func,nprocs,chunks,args=[],common_args=[]):
    """
    runs function in parallel for the list of args given in args
    """
    a_condition=all([len(a)==chunks for a in args if a not in common_args]) 
    if a_condition and len(args)>0:
        arranged_args=[]
        for pid in range(chunks):
            arglist =[]
            for a in args:
                if a not in common_args:
                    arglist.append(a[pid])
                else:
                    arglist.append(a)
            arranged_args.append(arglist)
        pool=mp.Pool(nprocs)
        results=[pool.apply_async(func=parallel_func, args=[func,pid,a]) for pid,a in zip(list(range(chunks)),arranged_args)]
        pool.close()
        pool.join()
        results=[r.get() for r in results]
        results.sort()
        results=[i[1] for i in results]
        return results
    else:
        pool=mp.Pool(nprocs)
        results=[pool.apply_async(func=parallel_func, args=[func,pid,args]) for pid in range(chunks)]
        pool.close()
        pool.join()
        results=[r.get() for r in results]
        results.sort()
        results=[i[1] for i in results]
        return results

def test_x_split_merge():
    """
    test for x_split and x_merge
    """
    x= np.random.rand(5,5)
    y= x_merge(x_split(x,2))
    assert np.all(x==y)
    x= np.random.rand(5,5,5)
    y= x_merge(x_split(x,2))
    assert np.all(x==y)

def test_y_split_merge():
    """
    test for y_split and y_merge
    """
    x= np.random.rand(5,5)
    y= y_merge(y_split(x,2))
    assert np.all(x==y)
    x= np.random.rand(5,5,5)
    y= y_merge(y_split(x,2))
    assert np.all(x==y)

def test_z_split_merge():
    """
    test for z_split and z_merge
    """
    x= np.random.rand(5,5,5)
    y= z_merge(z_split(x,2))
    assert np.all(x==y)


def add_array(x,y):
    """
    adds array used for testing
    """
    return x  + y
 
def list_operate(x):
    """
    opreates on list and returns its cube used for testing 
    """
    out=[]
    for i in x:
        time.sleep(0.1)
        out.append(i**3)
    return out

def test_list_parallel():
    """
    tests parallel implementation for list
    """    
    x=list(range(0,100))
    t0 = time.time()
    out_seq=list_operate(x)
    print("time taken for sequential simulation",time.time()-t0,'s')
    t0 = time.time()
    xx=list_split(x,4)
    out_parallel=run_parallel(list_operate,4,4,args=[xx])     
    out_parallel=list_merge(out_parallel)
    print("time taken for parallel simulation",time.time()-t0,'s')
    assert(out_seq==out_parallel)
    
def test_ndarray_parallel():
    """
    tests parallel ndarray implementation
    """
    x=np.ones((60,2500))
    y=2.*np.ones((60,2500))
    t0=time.time()
    zs=add_array(x,y)
    print("time taken for sequential simulation",time.time()-t0,'s')
    x=x_split(x,4)
    y=x_split(y,4)
    t0=time.time()
    results=run_parallel(add_array,4,4,args=[x,y])
    zp=x_merge(results)
    print("time taken for parallel simulation",time.time()-t0,'s')    
    assert np.all(zs==zp)

def call_imethod(method,instance,*args):
    return  getattr(instance,method)(*args)

class PrintX():
    def __init__(self,x):
        self.x = x 
    def print_x(self,prefix=''):
        time.sleep(0.5)
        print (prefix,self.x)
        
def test_call_imethod():
    ilist = [PrintX(1),PrintX(2),PrintX(3),PrintX(4),PrintX(5),PrintX(6),PrintX(7),PrintX(8)]
    t0=time.time()
    for i in ilist:
        i.print_x('hello')
    print("time taken for sequential simulation",time.time()-t0,'s')

    t0=time.time()
    run_parallel(call_imethod,8,8,
                    args=['print_x',ilist,'hello'],common_args=['print_x','hello'])
    print("time taken for parallel simulation",time.time()-t0,'s')    
                     
    
if __name__ == '__main__':
    test_x_split_merge()
    test_y_split_merge()
    test_z_split_merge()
    test_list_parallel()
    test_ndarray_parallel()
    test_call_imethod()
