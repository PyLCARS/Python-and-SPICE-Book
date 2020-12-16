#!/usr/bin/env python
# coding: utf-8

# This is not part of the book, this is a editing tool that is used to create markdownfiles for the jupyter noteboks and then convert them to Microsoft Word in order to use Word and Grammerly to do editing 

# In[1]:


from pathlib import Path


# In[2]:


path=Path('.')
notebooks=[]
for p in path.rglob("*.ipynb"):
    if '_build' in str(p.absolute()):
        continue
    
    elif 'ipynb_checkpoints' in str(p.absolute()):
        continue
    else:
        notebooks.append(p)
notebooks


# In[3]:


notebooks[2].relative_to(notebooks[2].parent.parent)


# In[4]:


for nb in notebooks:
    #base=str(nb.absolute())
    ipynb='./'+str(nb.relative_to(nb.parent.parent))
    base=ipynb[:-6]
    html=base+'.html'
    docx=base+'.docx'
    print(ipynb)
    
    #uncommment to use
    #!jupyter nbconvert --to html {ipynb}
    #!pandoc -s {html} -o {docx}
    
    


# In[ ]:




