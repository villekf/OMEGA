import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def volume3Dviewer(
    volume1: np.ndarray,
    scale = None,
    rotation = None,
    volume2 = None,
    useColorbar = False
):
    """
    Display an interactive Matplotlib slider to scroll through the third dimension of a 3D volume.

    Parameters
    ----------
    volume1 (np.ndarray):
        A 3D array of shape (X, Y, Z) or similar.
    scale : None, "fit", scalar, or (low, high)
        - None (default): normalize each slice independently (vmin/vmax per slice).
        - "fit": normalize all slices to the global min/max of volume1.
        - scalar: clip values below `scalar`, then normalize that clipped slice.
        - (low, high): clip values below low and above high, then normalize.
    rotation : None or [a, b, c]
        Which axis to slice along. Default None → [0,0,1] (Z-axis).  
        Other valid: [1,0,0] (X-axis) or [0,1,0] (Y-axis).
    volume2 : None or np.ndarray
        If provided, must be same shape as volume1. Overlays volume2 under its own colormap
        with a vertical slider controlling blend (alpha).
    """
    
    # Rotation
    if rotation not in [None, [0, 0, 1], [0, 1, 0], [1, 0, 0]]:
        raise ValueError("rotation must be one of None, [0,0,1], [1,0,0], [0,1,0]")
    elif rotation is None:
        axis = 2
    else:
        axis = rotation.index(1)
        
    # Second volume
    cmap1 = "gray"
    if volume2 is not None:
        if volume2.shape != volume1.shape:
            raise ValueError("volume2 must have the same shape as volume1")
        overlay = True
        cmap1 = "jet"
        cmap2 = "gray"
    else:
        overlay = False
    
    n_slices = volume1.shape[axis] # Number of slices in selected direction
    idx0 = n_slices // 2 # Default slice index
    
    # Set up figure and axes
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9, top=0.95)
    fig.canvas.manager.set_window_title('volume3Dviewer')

    # Helper to get the i‑th slice for chosen axis
    def get_slice_scale(vol, i):
        if axis == 0:
            s1 = vol[i, :, :]
        elif axis == 1:
            s1 = vol[:, i, :]
        else:
            s1 = vol[:, :, i]
        if scale is None:
            vmin, vmax = np.nanmin(s1), np.nanmax(s1)
        elif scale == "fit":
            vmin, vmax = np.nanmin(vol), np.nanmax(vol)
        elif np.isscalar(scale):
            vmin, vmax = scale, np.nanmax(s1)
        else:
            vmin, vmax = scale
        return s1, vmin, vmax

    def draw_slice(vol, idx, ax, cmap):
        sliceToDraw, vmin, vmax = get_slice_scale(vol, idx)
        img = ax.imshow(sliceToDraw, cmap=cmap, vmin=vmin, vmax=vmax)
        if overlay:
            img.set_alpha(0.5)
        return img
        
    # Draw initial slice
    img1 = draw_slice(volume1, idx0, ax, cmap1)
    if overlay:
        img2 = draw_slice(volume2, idx0, ax, cmap2)

    ax.axis('off')
    ax.set_title(f"Slice {idx0+1}/{n_slices}")

    # Horizontal slider for slice index
    ax_slice = plt.axes([.05, .05, .9, .05])
    slicer = Slider(ax_slice, '', 0, n_slices - 1, valinit=idx0, valstep=1)
    slicer.valtext.set_visible(False)
    # Vertical slider for blending
    if overlay:
        ax_alpha = plt.axes([.05, .15, .05, .8])
        alphar = Slider(ax_alpha, '', 0.0, 1.0, valinit=0.5, orientation='vertical')
        alphar.valtext.set_visible(False)
    
    if useColorbar:
        # Colorbar axes
        ax_cbar1 = plt.axes([.85, .15, .02, .8])
        _, vmin1, vmax1 = get_slice_scale(volume1, idx0)
        cbar1 = plt.colorbar(img1, cax=ax_cbar1, ticks=[vmin1, vmax1], ticklocation='left')
        #cbar1.ax.yaxis.set_ticks_position('left')
        if overlay:
            ax_cbar2 = plt.axes([.9, .15, .02, .8])
            _, vmin2, vmax2 = get_slice_scale(volume2, idx0)
            cbar2 = plt.colorbar(img2, cax=ax_cbar2, ticks=[vmin2, vmax2])
        
    # Update callback
    def _update(val):
        i = int(slicer.val)
        # get and scale slice1
        s1, vmin1, vmax1 = get_slice_scale(volume1, i)
        img1.set_data(s1)
        img1.set_clim(vmin1, vmax1)
        if not (np.isnan(vmin1) or np.isnan(vmax1)) and useColorbar:
            cbar1.set_ticks([vmin1, vmax1])

        if overlay:
            s2, v2min, v2max = get_slice_scale(volume2, i)
            img2.set_data(s2)
            img2.set_clim(v2min, v2max)
            if not (np.isnan(v2min) or np.isnan(v2max)) and useColorbar:
                cbar2.set_ticks([v2min, v2max])
            a = alphar.val
            img1.set_alpha(a)
            img2.set_alpha(1 - a)

        ax.set_title(f"Slice {i+1}/{n_slices}")
        fig.canvas.draw_idle()

    slicer.on_changed(_update)
    if overlay:
        alphar.on_changed(_update)

    plt.show()