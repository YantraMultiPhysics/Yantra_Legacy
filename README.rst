.. _PHREEQC: http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
.. _SCK-CEN: http://www.sckcen.be
.. _UGent: http://www.ugent.be/ea/structural-engineering/en/research/magnel
.. _LES: https://www.psi.ch/les/
.. _PSI: https://www.psi.ch/
.. _mailing: https://groups.google.com/forum/#!forum/yantra-mailing-list
.. _group: https://groups.google.com/forum/#!forum/yantra-users-queries
.. _TDM-MINGW: http://tdm-gcc.tdragon.net/
.. _Anaconda: https://www.anaconda.com/download/
.. _IPhreeqcPy: https://raviapatel.bitbucket.io/IPhreeqcPy/index.html

Yantra: A lattice Boltzmann Method based tool for multiscale/multiphysics simulations
=====================================================================================

*Yantra* is a lattice Boltzmann method based tool preliminary designed to perform pore-scale and multilevel pore-scale simulations. 
Lattice Boltzmann method derived from measoscopic Boltzmann equations provides certain advantages compared to traditional numerical 
methods when it comes to simulating complex processes at pore-scale. In *Yantra* additionally the transport module has been tightly
coupled with external geochemical solver `PHREEQC`_. Yantra also contains a special multilevel module which can handle simulations 
through porous media where due to lack of resolution both continuum and pore-scales coexist. 

Developer
++++++++++

Ravi A. Patel

email (personal): ravee.a.patel@gmail.com


Installation
++++++++++++
Currently yantra is in developmental stage. Therefore it advised not to install it
completely. Use following commands to build the FORTRAN extensions in place so that Yantra
can be ran from downloaded directory. At present python 3.12 and later versions are not supported.

Windows(may not work anymore)
-----------------------------

- Install 64 bit python distribution. Use of Anaconda is recommended.

- Install a FORTRAN compiler. Use of  TDM-MINGW distribution has been successfully tested.

- Install IPhreeqcPy. Follow its installation instructions.

- Download Yantra version 1.0.0-dev. Unzip the folder and run following commands  in TDM-GCC console  
	
	*python setup.py build_ext --inplace*

Linux (UBUNTU)
--------------

- Make sure  a Fortran compiler is available and numpy, matplotlib and scipy are installed 

- Install  requirements using requirements.txt file. Run following command in terminal

   *pip install -r requirements.txt*

- Download Yantra version 1.0.0-dev. Unzip the folder and run following command in terminal
	
	*python setup.py build_ext --inplace*

Detailed Documentation
+++++++++++++++++++++++
Coming soon

Support and Information
+++++++++++++++++++++++
In order to recieve news on developemental status of *Yantra* join `mailing`_ list.
For mailing list make sure you mark options allow email on every update to recieve updates

Please post your questions or queries regarding *Yantra* on this google `group`_

Found a bug in *Yantra*? report your bugs or suggestion use bitbucket issue tracker

Development history
+++++++++++++++++++
Development of *Yantra* was initated by the author during his PhD project 
which was carried out in joint collobration with Ghent university (`UGent`_)
and  Belgian Nuclear Reasearch Center (`SCK-CEN`_). At the end of PhD thesis (November 2015)
of the author *Yantra* version 0.1 was released.  Author continued to further develop *Yantra* during
November 2015 to July 2017 as part of his research at  Belgian Nuclear Reasearch Center (`SCK-CEN`_).
At the end of employment at `SCK-CEN`_ *Yantra* version 0.5 was released. From august the author started working at `LES`_ lab in PSI where he continues to develop and 
provide support for yantra. At present *Yantra* is being redesigned and a new version 1.0.0 is due to be released
in January 2018

License and Terms of use
++++++++++++++++++++++++

*Yantra* is a free software: you can redistribute it and/or modify it 
under the terms of the GNU  General Public License as published by the
Free Software Foundation, version 3 of the License. This program is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU Lesser General Public License for more details. You should have 
received a copy of the GNU Lesser General Public License along with this program.
If not, see `<http://www.gnu.org/licenses/>`_.

Institutes
++++++++++

.. image:: /images/institutes_involved.jpg
   :height: 100px
   :width: 200 px
   :scale: 50 %
   :align: center
   
