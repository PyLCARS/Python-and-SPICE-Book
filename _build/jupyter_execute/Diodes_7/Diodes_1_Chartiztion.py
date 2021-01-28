#!/usr/bin/env python
# coding: utf-8

# NOT STARTED

# In[1]:


from PySpice.Spice.Library import SpiceLibrary as spl


# In[2]:


lib=spl('KiCad-Spice-Library/Models')


# In[3]:


len(lib._models.keys())


# In[4]:


lib.search('.*2222.*')


# In[ ]:




