#!/usr/bin/env python
# coding: utf-8

# 
# # Dependency install instruction
# Disclaimer this was all done with Ubuntu 20.04.1 LTS via the WSL2 for Windows 20.
# The ngspice version is ngspice-32 that was installed via
# ```
# Sudo apt-get install ngspice
# ```
# The SKiDl version that is being used is the authors personal fork of SKiDl 0.0.31.dev0 with some modifications that are being added to the in development branch of SKiDl as of 20201214 and can be got via
# ```
# pip install git+https://github.com/GProtoZeroW/skidl.git
# ```
# When the in development branch of version of SKiDl become the master one I will update this and the book to reflect the said SKiDl release. FYI the SKiDl main repo is located at https://github.com/xesscorp/skidl
# 
# The PySpice version that is being used is the authors personal fork of PySpice 1.4.3 with some modifications that have yet to be pushed to the main PySpice repo which is found at https://github.com/FabriceSalvaire/PySpice . When the changes are incorporated into the main version repo of PySpice I will update this and the book to reflect the said PySpice release enhancements for this book.  You can install the fork with the needed modifications for this book to work via 
# ```
# pip install git+https://github.com/GProtoZeroW/PySpice.git
# ```
# The version of Lcapy being used is 0.70 and comes strait from it source repo https://github.com/mph-/lcapy and can be install via 
# ```
# pip install lcapy
# ```
# All other python libraries being used are part of the standard python scientific stack as part of the Anaconda distribution
# 

# In[1]:


#tool to log notebook internals
#https://github.com/jrjohansson/version_information
get_ipython().run_line_magic('load_ext', 'version_information')
get_ipython().run_line_magic('version_information', 'skidl, PySpice,lcapy, sympy, numpy, matplotlib, pandas, scipy')


# # Book Specific Code
# The dependent code developed in writing this book can be found inside each chapters folder containing alongside the source Jupyter notebooks to build the chapter in a python file named:
# ```
# <chapter_name>_<chapter_number>_Codes.py
# ```
# Where each entry in the .py file will state what chapter number and what section number of the jupyter notebook where the code came from
# 

# 
# ```{toctree}
# :hidden:
# :titlesonly:
# 
# 
# DC_1/DC_1
# ```
# 
