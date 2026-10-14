"""Pack full-FOV images into the projector's shifted resolution volumes."""
import numpy as np


def _resize_linear(image, size):
    """Pixel-center linear resize with antialiasing and symmetric boundaries."""
    from scipy.sparse import csr_matrix

    result = image.astype(np.float32)
    for axis, output_size in enumerate(size):
        input_size = result.shape[axis]
        if input_size == output_size:
            continue
        scale = float(output_size) / input_size
        kernel_scale = min(scale, 1.)
        centers = (np.arange(output_size) + 0.5) / scale - 0.5
        left = np.floor(centers - 1. / kernel_scale).astype(np.int64)
        indices = left[:, None] + np.arange(int(np.ceil(2. / kernel_scale)) + 2)
        weights = np.maximum(0., 1. - abs(indices - centers[:, None]) * kernel_scale)
        weights /= weights.sum(axis=1, keepdims=True)
        reflected = indices % (2 * input_size)
        reflected = np.where(reflected < input_size, reflected, 2 * input_size - 1 - reflected)
        rows = np.broadcast_to(np.arange(output_size)[:, None], indices.shape)
        matrix = csr_matrix((weights.ravel().astype(np.float32), (rows.ravel(), reflected.ravel())),
                            shape=(output_size, input_size))
        moved = np.moveaxis(result, axis, 0)
        resized = (matrix @ moved.reshape(input_size, -1)).reshape((output_size, *moved.shape[1:]))
        result = np.moveaxis(resized, 0, axis)
    return result


def pack_multiresolution(image, options, mask=False):
    from scipy.ndimage import zoom

    sizes = np.column_stack((options.Nx, options.Ny, options.Nz)).astype(np.int64)
    scale = float(options.multiResolutionScale)
    main = np.floor(sizes[0] * scale + 0.5).astype(np.int64)
    count = len(sizes)
    if count == 7:
        before = np.array([sizes[3, 0], sizes[5, 1], sizes[1, 2]])
        after = np.array([sizes[4, 0], sizes[6, 1], sizes[2, 2]])
        starts = [(before[0], before[1], 0), (before[0], before[1], before[2]+main[2]),
                  (0, 0, 0), (before[0]+main[0], 0, 0),
                  (before[0], 0, 0), (before[0], before[1]+main[1], 0)]
    elif count == 5:
        before = np.array([sizes[1, 0], sizes[3, 1], 0])
        after = np.array([sizes[2, 0], sizes[4, 1], 0])
        starts = [(0, 0, 0), (before[0]+main[0], 0, 0),
                  (before[0], 0, 0), (before[0], before[1]+main[1], 0)]
    elif count == 3:
        before = np.array([0, 0, sizes[1, 2]])
        after = np.array([0, 0, sizes[2, 2]])
        starts = [(0, 0, 0), (0, 0, before[2]+main[2])]
    else:
        raise ValueError('Multi-resolution packing requires 3, 5, or 7 volumes.')
    low_size = np.maximum(before + main + after, sizes[1:].max(axis=0))
    high_start = np.floor(before / scale + 0.5).astype(np.int64)
    high_size = high_start + sizes[0] + np.floor(after / scale + 0.5).astype(np.int64)
    image = np.asarray(image)
    if image.ndim == 1:
        image = image.reshape((int(options.NxFull), int(options.NyFull), int(options.NzFull)), order='F')
    if mask and image.ndim == 2:
        image = np.repeat(image[:, :, None], int(options.NzFull), axis=2)
    if image.ndim != 3:
        raise ValueError('Multi-resolution input must be a full-FOV 3D image.')
    def resized(size):
        if mask:
            return zoom(image.astype(np.float32), size / np.asarray(image.shape),
                        order=0, mode='nearest', grid_mode=True, prefilter=False)
        return _resize_linear(image, size)
    def crop(volume, start, size):
        result = volume[tuple(slice(int(a), int(a+n)) for a, n in zip(start, size))]
        if tuple(result.shape) != tuple(size):
            raise ValueError('Multi-resolution crop does not match its volume dimensions.')
        return result.ravel(order='F')
    high, low = resized(high_size), resized(low_size)
    volumes = [crop(high, high_start, sizes[0])]
    volumes.extend(crop(low, start, size) for start, size in zip(starts, sizes[1:]))
    packed = np.concatenate(volumes)
    return (packed > 0).astype(np.uint8) if mask else packed.astype(np.float32)
