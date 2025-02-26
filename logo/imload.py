# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import imageio
media=imageio.imread('logo_media.png')
media = 1-(media [:,:,2]>0)
media[media==0]=-1
import matplotlib.pyplot as plt
plt.imshow(media)
plt.show()