# -*- coding: utf-8 -*-
import imageio
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM

images = []
for i in range(5000,250000,5000):
    drawing = svg2rlg('streamline_itr_%s.svg'%i)
    renderPM.drawToFile(drawing, "file.png", fmt="PNG")
    images.append(imageio.imread('file.png'))
imageio.mimsave('movie1.gif', images)
