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

from setuptools import Extension,setup
from setuptools.command.build_ext import build_ext
from os import listdir
from os.path import join, isfile, splitext
import sys,os
import shutil
class f2py_extensions(Extension):
    def __init__(self, name, sysargs, source, module, module_name,module_loc,fname):
        Extension.__init__(self, name = name, sources = [])
        self.dirs = source
        self.module = module
        self.module_name = module_name
        self.module_loc = module_loc
        self.fnames = fname
        self.sysargs = sysargs
        #compiler 
        v=[v for v in sysargs if v.startswith('--fcompiler')]
        if len(v)>0:
            self.compiler = v[0].split("=")[1].lower()
        else:
            self.compiler = "gnu95"
        self.complier_args, self.link_args = self.get_compiler_options(self.compiler)
    
    @staticmethod
    def get_compiler_options(compiler:str)->tuple:  
        """
        Get compiler options based on the compiler name
        """
        compiler_options = {
        "gnu": ["g95", "gnu95", "mingw32"],
        "intel": ["intel", "intelm", "intelem", "intelv", "intelvem"],
        "pg": ["pg", "pgfortran", "pgf90", "pgf95"]
        }

        compile_args = {
            "gnu": ["-fopenmp", "-fPIC", "-O3", "-fbounds-check", "-mtune=native"],
            "intel": ["-qopenmp", "-O3", "-funroll-loops"],
            "pg": ["-mp"],
            "default": ["-fopenmp", "-O3"]
        }

        link_args = {
            "gnu": ["-lgomp"],
            "intel": ["-liomp5"],
            "pg": [],
            "default": ["-lgomp"]
        }
        
        for key, compilers in compiler_options.items():
            if compiler in compilers:
                return compile_args[key], link_args[key]
        return compile_args["default"], link_args["default"]


class f2py_build_ext(build_ext):
    """
    builts f2py extension modules
    """
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self,ext:f2py_extensions):
        compiler_args = " ".join(ext.complier_args)
        link_args = " ".join(ext.link_args)
        compiler_name = ext.compiler
        for mloc, mn, fn in zip(ext.module_loc, ext.module_name, ext.fnames):
            print(f"Building {mn} in {mloc}")
            os.system(
            f"cd {mloc};f2py -c {fn} -m {mn} --fcompiler={compiler_name} --f90flags='{compiler_args}' {link_args}"
            )


def list_sources_and_modules(parent_folder,parent_package,extensions = ['py'],
        prefix = '',exclude = [],result = None) -> dict:
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
        result = {'source': [],'module':[], 'module_name': [], 'module_loc': [], "fname": []}
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
                    result['module_name'].append(name)
                    result['module_loc'].append(os.path.abspath(parent_folder))
                    result['source'].append(join(parent_folder,fname))
                    result['fname'].append(fname)
            else:
                result = list_sources_and_modules(path,
                                                  parent_package+'.'+fname,
                                                  extensions, prefix, exclude,
                                                  result)
    return result

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
    ext_modules = f2py_extensions("fortran_modules", args,**fort)
    #setup
    setup(
        name = 'yantra',
        version=version,
        author = 'Ravi A. Patel',
        author_email = 'ravee.a.patel@gmail.com',
        py_modules= py['module'],
        ext_modules = [ext_modules],
        cmdclass={'build_ext': f2py_build_ext},
        license='GPL V3 and later',
         )
    print ('+++Setup complete.')    
if __name__ == '__main__':
    #remove existing build directory
    if os.path.isdir('build'):
        shutil.rmtree('build')
    #runs setup
    run_setup(sys.argv)


