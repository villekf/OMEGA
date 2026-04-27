import numpy as np

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
    
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    
    # Rotation
    if rotation not in [None, [0, 0, 1], [0, 1, 0], [1, 0, 0]]:
        raise ValueError("rotation must be one of None, [0,0,1], [1,0,0], [0,1,0]")
    elif rotation is None:
        axis = 2
    else:
        axis = rotation.index(1)

    # Apply rotation permutation to 3D/4D
    if rotation is not None:
        if axis == 0:
            if volume1.ndim == 4:
                volume1 = np.flip(np.transpose(volume1, (0, 2, 1, 3)), axis=0)
            else:
                volume1 = np.flip(np.transpose(volume1, (0, 2, 1)), axis=0)
        elif axis == 1:
            if volume1.ndim == 4:
                volume1 = np.flip(np.transpose(volume1, (2, 1, 0, 3)), axis=0)
            else:
                volume1 = np.flip(np.transpose(volume1, (2, 1, 0)), axis=0)

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

    # 4D check
    has_4d = volume1.ndim == 4 and volume1.shape[3] > 1
    n_slices = volume1.shape[2]
    n_dim4 = volume1.shape[3] if has_4d else 1
    idx0 = n_slices // 2
    idx4_0 = 0

    # Set up figure and axes
    fig, ax = plt.subplots()
    # Leave more bottom space when 4D (two sliders)
    bottom_margin = 0.18 if has_4d else 0.12
    plt.subplots_adjust(left=0.1, bottom=bottom_margin, right=0.9, top=0.95)
    fig.canvas.manager.set_window_title('volume3Dviewer')

    def get_slice(vol, i, i4):
        if has_4d:
            s = vol[:, :, i, i4]
        else:
            s = vol[:, :, i]
        if scale is None:
            vmin, vmax = np.nanmin(s), np.nanmax(s)
        elif scale == "fit":
            vmin, vmax = np.nanmin(vol), np.nanmax(vol)
        elif np.isscalar(scale):
            vmin, vmax = scale, np.nanmax(s)
        else:
            vmin, vmax = scale
        return s, vmin, vmax

    def draw_slice(vol, idx, idx4, ax, cmap):
        sliceToDraw, vmin, vmax = get_slice(vol, idx, idx4)
        img = ax.imshow(sliceToDraw, cmap=cmap, vmin=vmin, vmax=vmax)
        if overlay:
            img.set_alpha(0.5)
        return img

    # Draw initial slice
    img1 = draw_slice(volume1, idx0, idx4_0, ax, cmap1)
    if overlay:
        img2 = draw_slice(volume2, idx0, idx4_0, ax, cmap2)

    ax.axis('off')

    def make_title(i, i4):
        if has_4d:
            return f"Slice {i+1}/{n_slices}  (T: {i4+1}/{n_dim4})"
        return f"Slice {i+1}/{n_slices}"

    ax.set_title(make_title(idx0, idx4_0))

    # Horizontal slider for slice index (Z)
    # Position: if 4D, shift up to make room for the second slider
    slice_slider_bottom = 0.09 if has_4d else 0.04
    ax_slice = plt.axes([0.05, slice_slider_bottom, 0.9, 0.04])
    slicer = Slider(ax_slice, 'Z', 0, n_slices - 1, valinit=idx0, valstep=1)
    slicer.valtext.set_visible(False)

    # Horizontal slider for 4th dimension (T)
    if has_4d:
        ax_dim4 = plt.axes([0.05, 0.03, 0.9, 0.04])
        slider_dim4 = Slider(ax_dim4, 'T', 0, n_dim4 - 1, valinit=idx4_0, valstep=1)
        slider_dim4.valtext.set_visible(False)

    # Vertical slider for blending (overlay)
    if overlay:
        ax_alpha = plt.axes([0.05, 0.20, 0.05, 0.75])
        alphar = Slider(ax_alpha, '', 0.0, 1.0, valinit=0.5, orientation='vertical')
        alphar.valtext.set_visible(False)

    if useColorbar:
        ax_cbar1 = plt.axes([0.85, 0.20, 0.02, 0.75])
        _, vmin1, vmax1 = get_slice(volume1, idx0, idx4_0)
        cbar1 = plt.colorbar(img1, cax=ax_cbar1, ticks=[vmin1, vmax1], ticklocation='left')
        if overlay:
            ax_cbar2 = plt.axes([0.9, 0.20, 0.02, 0.75])
            _, vmin2, vmax2 = get_slice(volume2, idx0, idx4_0)
            cbar2 = plt.colorbar(img2, cax=ax_cbar2, ticks=[vmin2, vmax2])

    def _update(val):
        i = int(slicer.val)
        i4 = int(slider_dim4.val) if has_4d else 0

        s1, vmin1, vmax1 = get_slice(volume1, i, i4)
        img1.set_data(s1)
        img1.set_clim(vmin1, vmax1)
        if useColorbar and not (np.isnan(vmin1) or np.isnan(vmax1)):
            cbar1.set_ticks([vmin1, vmax1])

        if overlay:
            s2, v2min, v2max = get_slice(volume2, i, i4)
            img2.set_data(s2)
            img2.set_clim(v2min, v2max)
            if useColorbar and not (np.isnan(v2min) or np.isnan(v2max)):
                cbar2.set_ticks([v2min, v2max])
            a = alphar.val
            img1.set_alpha(a)
            img2.set_alpha(1 - a)

        ax.set_title(make_title(i, i4))
        fig.canvas.draw_idle()

    slicer.on_changed(_update)
    if has_4d:
        slider_dim4.on_changed(_update)
    if overlay:
        alphar.on_changed(_update)

    # Prevent garbage collection of widget objects in interactive environments
    fig._widgets = [slicer]
    if has_4d:
        fig._widgets.append(slider_dim4)
    if overlay:
        fig._widgets.append(alphar)

    plt.show()