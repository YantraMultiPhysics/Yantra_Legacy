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
#PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#=======================================================================================

from numpy.distutils.core import Extension
from numpy.distutils.core import setup
from os import listdir
from os.path import join, isfile, splitext
import sys,os
import shutil
def list_sources_and_modules(parent_folder,parent_package,extensions = ['py'],
        prefix = '',exclude = [],result = None):
    """
    creates list of sources and module names by going through the parent folder
    
    Parameters
    ----------
    parent_folder: str
        name of the parent folder to look into
    parent_package: str
        name of the parent package 
    extensions: list
        list of file extensions to look for
    prefix: str
        prefix if any to be included in the name
    exclude: str
        path to be excluded for the search
    result: dict, Default None
        dict of the existing output of this function
        
    returns
    -------
        dict
            dictionary containing keys 'module' and 'source' corresponding to the list of modules name and their sources respectively
    """
    if result == None:
        result = {'source': [], 'module': []}
    for fname in listdir(parent_folder):
        path = join(parent_folder,fname)
        acceptable = path not in exclude and fname not in exclude
        acceptable = acceptable and "_ignore" not in fname
        if acceptable:
            if isfile(path):
                name, ext = splitext(fname)
                if ext in extensions:
                    name = prefix + name
                    result['module'].append(parent_package+'.'+name)
                    result['source'].append(join(parent_folder,fname))
            else:
                result = list_sources_and_modules(path,
                                                  parent_package+'.'+fname,
                                                  extensions, prefix, exclude,
                                                  result)
    return result

def get_ext_modules(name,source,args = []):
    """
    Gets list of extension instances from module name and its corresponding fortran source file
    
    Parameters
    ----------
    name: list
        list of module names
    source: list
        list of path of the source files
    args: list
        
    Returns
    -------
    list
        list of extension instances
    """
    v=[v for v in args if v.startswith('--fcompiler')]
    if len(v)>0:
        compiler = v[0].split("=")[1].lower()
        print(compiler)
        gnu_compiler = ["g95",'gnu95',"mingw32"]
        intel_compiler = ["intel","intelm","intelem","intelv","intelvem"]
        pg_compiler=["pg","pgfortran", "pgf90", "pgf95"]
        if compiler in gnu_compiler:
            extra_f90_compile_args = ["-fopenmp", "-fPIC", "-O3", "-fbounds-check",
                                      "-mtune=native"]
#            extra_f90_compile_args = ["-fopenmp", "-O3"]
            extra_link_args = ["-lgomp"]            
        elif compiler in intel_compiler:
#            extra_f90_compile_args = ["-openmp", "-fPIC", "-xHost", "-O3", "-ipo",
#                                      "-funroll-loops", "-heap-arrays", "-mcmodel=medium"]
            extra_f90_compile_args = ["-qopenmp", "-O3", "-funroll-loops"]
            extra_link_args = ["-liomp5"]
        elif compiler in pg_compiler:
            extra_f90_compile_args = ["-mp"]
            extra_link_args=[]
        else:
            raise ValueError("Support for this compiler not included currently in setup file contact the developer")
    else:
        extra_f90_compile_args = ["-fopenmp","-O3"]
        extra_link_args = ["-lgomp"]           
    ext_modules = []
    for (module, src) in zip(name,source):
        fext = Extension(name = module,
                         sources = [src],
                         extra_f90_compile_args = extra_f90_compile_args,
                         extra_link_args = extra_link_args,
                         )
        ext_modules.append(fext)
    return ext_modules
    
def run_setup(args):
    """
    runs setup to install yantra
    """
    #write version
    major = 1
    minor = 0
    patch = 0
    release = '-dev'
    version = '{0}.{1}.{2}{3}'.format(major,minor,patch,release)
    fname = os.path.join('yantra','version.py')
    with open(fname, 'w') as f:
        f.write("version = '{0}'".format(version))
    #get list of fortran and python modules
    fort = list_sources_and_modules('yantra', 'yantra', ['.f90'],
                                    prefix = '', exclude = [])
    py = list_sources_and_modules('yantra', 'yantra', ['.py'],
                                  prefix = '', exclude = [])
    #build extension instances for fortran modules
    ext_modules = get_ext_modules(fort['module'],fort['source'],args)
    #setup
    setup(
        name = 'yantra',
        version=version,
        author = 'Ravi A. Patel',
        author_email = 'ravee.a.patel@gmail.com',
        py_modules= py['module'],
        ext_modules = ext_modules,
        license='GPL V3 and later',
         )
    print ('+++Setup complete.')    
if __name__ == '__main__':
    #remove existing build directory
    if os.path.isdir('build'):
        shutil.rmtree('build')
    #runs setup
    run_setup(sys.argv)


