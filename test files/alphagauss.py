# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 12:23:15 2019

@author: wsraymon
"""

# sphinx_gallery_thumbnail_number = 3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def normal_pdf(x, mean, var):
    return np.exp(-(x - mean)**2 / (2*var))


# Generate the space in which the blobs will live
xmin, xmax, ymin, ymax = (0, 100, 0, 100)
n_bins = 100
xx = np.linspace(xmin, xmax, n_bins)
yy = np.linspace(ymin, ymax, n_bins)

# Generate the blobs. The range of the values is roughly -.0002 to .0002
means_high = [20, 50]
means_low = [50, 60]
var = [150, 200]

gauss_x_high = normal_pdf(xx, means_high[0], var[0])
gauss_y_high = normal_pdf(yy, means_high[1], var[0])

gauss_x_low = normal_pdf(xx, means_low[0], var[1])
gauss_y_low = normal_pdf(yy, means_low[1], var[1])

weights_high = np.array(np.meshgrid(gauss_x_high, gauss_y_high)).prod(0)
weights_low = -1 * np.array(np.meshgrid(gauss_x_low, gauss_y_low)).prod(0)
weights = weights_high + weights_low

# We'll also create a grey background into which the pixels will fade
greys = np.empty(weights.shape + (3,), dtype=np.uint8)
greys.fill(0)

# First we'll plot these blobs using only ``imshow``.
vmax = np.abs(weights).max()
vmin = -vmax
cmap = plt.cm.viridis

fig, ax = plt.subplots()
ax.imshow(greys)
ax.imshow(weights, extent=(xmin, xmax, ymin, ymax), cmap=cmap)
ax.set_axis_off()


# Create an alpha channel based on weight values
# Any value whose absolute value is > .0001 will have zero transparency
alphas = Normalize(0, 1, clip=True)(np.abs(weights))
#alphas = np.clip(alphas,0, 1)  # alpha value clipped at the bottom at .4

# Normalize the colors b/w 0 and 1, we'll then pass an MxNx4 array to imshow
colors = Normalize(vmin, vmax)(weights)
colors = cmap(colors)
print(colors.shape)

# Now set the alpha channel to the one we created above
colors[..., -1] = alphas
print(colors.shape)

# Create the figure and image
# Note that the absolute values may be slightly different
fig, ax = plt.subplots()
#ax.imshow(greys)
#ax.imshow(colors, extent=(xmin, xmax, ymin, ymax))

greys = np.empty(weights.shape + (3,), dtype=np.uint8)
greys.fill(0)
ax.imshow(greys)


n_bins = 100
xx = np.linspace(xmin, xmax, n_bins)
yy = np.linspace(ymin, ymax, n_bins)        


weights = np.array(np.meshgrid(xx, yy)).prod(0)
greys = np.empty(weights.shape + (4,), dtype=np.float32)
greys.fill(100)  
 

greys = greys + np.random.poisson(greys)
alphas = Normalize(0, 1, clip=True)(np.abs(greys))
alphas[:,:,:] = .5


greys = Normalize(0, 256)(greys)
greys[..., -1] = .3
ax.imshow(colors,extent=(xmin, xmax, ymin, ymax))
ax.imshow(greys,extent=(xmin, xmax, ymin, ymax))

# Add contour lines to further highlight different levels.
#ax.contour(weights[::-1], levels=[-.1, .1], colors='k', linestyles='-')
#ax.set_axis_off()
#plt.show()

#ax.contour(weights[::-1], levels=[-.0000001, .0000001], colors='k', linestyles='-')
ax.set_axis_off()
plt.show()









import time


n_bins = 5000
xx = np.linspace(xmin, xmax, n_bins)
yy = np.linspace(ymin, ymax, n_bins)        


weights = np.array(np.meshgrid(xx, yy)).prod(0)

def plot_gauss(centers,rs,weight):
    st = time.time()
    xmin, xmax, ymin, ymax = (-60, 60, -60, 60)
    n_bins = 100
    xx = np.linspace(xmin, xmax, n_bins)
    yy = np.linspace(ymin, ymax, n_bins)
    fig, ax = plt.subplots()
  
        
    i = 0
    for center in centers:
        weightsx = normal_pdf(xx, center[0], 10)
        weightsy = normal_pdf(yy, -center[1], 10)
        
        weights = np.array(np.meshgrid(weightsx, weightsy)).prod(0)
        if i == 0:
            greys = np.empty(weights.shape + (3,), dtype=np.uint8)
            greys.fill(0)  
            i+=1
            ax.imshow(greys)
    
        vmax = np.abs(weights).max()
        vmin = -vmax
        cmap = plt.cm.viridis
        
        alphas = Normalize(0, 1, clip=True)(np.abs(weights))
        
        colors = Normalize(vmin, vmax)(weights)
        colors = cmap(colors)
        
        # Now set the alpha channel to the one we created above
        colors[..., -1] = alphas
        
        # Create the figure and image
        # Note that the absolute values may be slightly different
        
       
        ax.imshow(colors, extent=(xmin, xmax, ymin, ymax))
        
    print(time.time()-st)
