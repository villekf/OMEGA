# -*- coding: utf-8 -*-
"""Native PyTorch-MPS bridge for OMEGA projector kernels."""

from __future__ import annotations

import hashlib
import os
import re
import shlex
import struct
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCALAR_KERNEL_PARAMS_SIZE = 352

_OFFSETS = {
    "nRowsD": 0,
    "nColsD": 4,
    "dPitch": 8,
    "dL": 16,
    "global_factor": 20,
    "epps": 24,
    "det_per_ring": 28,
    "sigma_x": 32,
    "coneOfResponseStdCoeffA": 36,
    "coneOfResponseStdCoeffB": 40,
    "coneOfResponseStdCoeffC": 44,
    "tube_width": 48,
    "cylRadiusProj3": 52,
    "bmin": 56,
    "bmax": 60,
    "Vmax": 64,
    "rings": 68,
    "helicalRadius": 72,
    "ellipseCenter": 80,
    "ellipseRadii": 96,
    "ellipsePower": 112,
    "d_N": 128,
    "b": 144,
    "dSize5": 160,
    "kerroin4": 168,
    "DSC": 172,
    "d": 176,
    "d_Scale4": 192,
    "d_Scale5": 208,
    "d_bmax": 224,
    "orthWidth": 240,
    "nProjections": 248,
    "no_norm": 256,
    "m_size": 264,
    "currentSubset": 272,
    "aa": 276,
    "N_PDHG": 288,
    "epps_PDHG": 304,
    "theta_PDHG": 308,
    "tau_PDHG": 312,
    "enforcePositivity_PDHG": 316,
    "N_rotate": 320,
    "cosa_rotate": 336,
    "sina_rotate": 340,
}


def _empty_value(value: Any) -> bool:
    try:
        return np.asarray(value).size == 0
    except Exception:
        return value is None


def _pack_scalar_kernel_params(self: Any, timestep: int, subset: int, volume: int) -> bytes:
    blob = bytearray(SCALAR_KERNEL_PARAMS_SIZE)

    def u32(name: str, value: int) -> None:
        struct.pack_into('<I', blob, _OFFSETS[name], int(value))

    def i32(name: str, value: int) -> None:
        struct.pack_into('<i', blob, _OFFSETS[name], int(value))

    def i64(name: str, value: int) -> None:
        struct.pack_into('<q', blob, _OFFSETS[name], int(value))

    def u64(name: str, value: int) -> None:
        struct.pack_into('<Q', blob, _OFFSETS[name], int(value))

    def u8(name: str, value: int) -> None:
        struct.pack_into('<B', blob, _OFFSETS[name], int(value))

    def f32(name: str, value: float) -> None:
        struct.pack_into('<f', blob, _OFFSETS[name], float(value))

    def f2(name: str, values: Iterable[float]) -> None:
        a, b = values
        struct.pack_into('<2f', blob, _OFFSETS[name], float(a), float(b))

    def f3(name: str, values: Iterable[float]) -> None:
        a, b, c = values
        struct.pack_into('<4f', blob, _OFFSETS[name], float(a), float(b), float(c), 0.0)

    def u3(name: str, values: Iterable[int]) -> None:
        a, b, c = values
        struct.pack_into('<4I', blob, _OFFSETS[name], int(a), int(b), int(c), 0)

    def i3(name: str, values: Iterable[int]) -> None:
        a, b, c = values
        struct.pack_into('<4i', blob, _OFFSETS[name], int(a), int(b), int(c), 0)

    nx = self.Nx[volume]
    ny = self.Ny[volume]
    nz = self.Nz[volume]
    dx = self.dx[volume]
    dy = self.dy[volume]
    dz = self.dz[volume]
    bx = self.bx[volume]
    by = self.by[volume]
    bz = self.bz[volume]

    u32('nRowsD', int(self.nRowsD))
    u32('nColsD', int(self.nColsD))
    f2('dPitch', (float(self.dPitchX), float(self.dPitchY)))
    f32('dL', float(self.dL))
    f32('global_factor', float(self.global_factor))
    f32('epps', float(self.epps))
    u32('det_per_ring', int(self.det_per_ring))
    f32('sigma_x', float(self.sigma_x))
    f32('coneOfResponseStdCoeffA', float(self.coneOfResponseStdCoeffA))
    f32('coneOfResponseStdCoeffB', float(self.coneOfResponseStdCoeffB))
    f32('coneOfResponseStdCoeffC', float(self.coneOfResponseStdCoeffC))
    f32('tube_width', float(self.tube_width_z))
    f32('cylRadiusProj3', float(self.tube_radius))
    f32('bmin', float(self.bmin))
    f32('bmax', float(self.bmax))
    f32('Vmax', float(self.Vmax))
    u32('rings', int(self.rings))
    f32('helicalRadius', float(self.helicalRadius))
    f3('ellipseCenter', (
        float(self.ellipseCenterX),
        float(self.ellipseCenterY),
        float(self.ellipseCenterZ),
    ))
    f3('ellipseRadii', (
        float(self.ellipseRadiusX),
        float(self.ellipseRadiusY),
        float(self.ellipseRadiusZ),
    ))
    f32('ellipsePower', float(self.ellipsePower))

    u3('d_N', (nx, ny, nz))
    f3('b', (bx, by, bz))
    if self.FPType == 5 or self.BPType == 5:
        f2('dSize5', (self.dSizeX[volume], self.dSizeY[volume]))
        f3('d_Scale5', (self.dScaleX[volume], self.dScaleY[volume], self.dScaleZ[volume]))
    if self.BPType == 4:
        f32('kerroin4', float(self.kerroin[volume]))
    f32('DSC', float(getattr(self, 'DSC', 0.0)))
    f3('d', (dx, dy, dz))
    if self.FPType == 4:
        f3('d_Scale4', (self.dScaleX4[volume], self.dScaleY4[volume], self.dScaleZ4[volume]))
    f3('d_bmax', (
        bx + nx * dx,
        by + ny * dy,
        bz + nz * dz,
    ))
    f32('orthWidth', float(self.tube_width_z))
    i64('nProjections', int(self.nProjSubset[timestep, subset]))
    u8('no_norm', int(self.no_norm))
    u64('m_size', int(self.nMeasSubset[timestep, subset]))
    u32('currentSubset', subset)
    i32('aa', volume)
    i3('N_PDHG', (0, 0, 0))
    f32('epps_PDHG', 0.0)
    f32('theta_PDHG', 0.0)
    f32('tau_PDHG', 0.0)
    u8('enforcePositivity_PDHG', 0)
    i3('N_rotate', (0, 0, 0))
    f32('cosa_rotate', 0.0)
    f32('sina_rotate', 0.0)
    return bytes(blob)


def _mps_tensor_from_numpy(torch: Any, value: Any, dtype: Any) -> Any:
    arr = np.asarray(value, dtype=dtype)
    return torch.as_tensor(np.ascontiguousarray(arr), device='mps')


def _validate_configuration(self: Any) -> None:
    errors: list[str] = []
    if self.FPType not in (1, 2, 3, 6):
        errors.append(f'forward projector type {self.FPType!r} is unsupported')
    if self.BPType not in (1, 2, 3, 4, 6):
        errors.append(f'backprojector type {self.BPType!r} is unsupported')
    if (self.FPType == 6 or self.BPType == 6) and not self.SPECT:
        errors.append('projector type 6 is only supported with SPECT data')
    if self.use_32bit_atomics or self.use_64bit_atomics:
        errors.append('integer accumulation is unsupported by the Metal/MPS custom-operator path')
    if self.use_psf:
        errors.append('PSF convolution is unsupported by the Metal/MPS custom-operator path')
    if errors:
        raise ValueError('Metal/MPS configuration: ' + '; '.join(errors))


_SHADER_CACHE: dict[str, Any] = {}
_LOCAL_INCLUDE_RE = re.compile(r'^\s*#\s*include\s*"([^"]+)"\s*(?://.*)?$')


def _iter_option_tokens(options: Iterable[Any]) -> Iterable[str]:
    for option in options:
        text = str(option).strip()
        if not text:
            continue
        try:
            tokens = shlex.split(text)
        except ValueError:
            tokens = text.split()
        yield from (token.strip() for token in tokens)


def _macro_preamble(options: Iterable[Any]) -> str:
    definitions: dict[str, str | None] = {}
    order: list[str] = []
    for token in _iter_option_tokens(options):
        if not token.startswith('-D') or len(token) <= 2:
            continue
        body = token[2:]
        if '=' in body:
            name, value = body.split('=', 1)
        else:
            name, value = body, None
        name = name.strip()
        if not name or name in {'PYTHON', 'USEIMAGES'}:
            continue
        if name not in definitions:
            order.append(name)
        definitions[name] = value.strip() if value is not None else None
    if 'METAL' not in definitions:
        order.insert(0, 'METAL')
        definitions['METAL'] = None
    lines = ['// Generated from OMEGA projector settings.']
    for name in order:
        value = definitions[name]
        lines.append(f'#define {name} 1' if value in (None, '') else f'#define {name} {value}')
    return '\n'.join(lines) + '\n\n'


def _resolve_local_include(name: str, current_dir: Path, search_dirs: tuple[Path, ...]) -> Path:
    candidates = [current_dir / name, *(directory / name for directory in search_dirs)]
    for candidate in candidates:
        candidate = candidate.resolve()
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f'Unable to resolve Metal include {name!r}; searched: {candidates}')


def _inline_local_includes(
    source: str,
    *,
    current_dir: Path,
    search_dirs: tuple[Path, ...],
    active_stack: tuple[Path, ...] = (),
) -> str:
    output: list[str] = []
    for line in source.splitlines():
        match = _LOCAL_INCLUDE_RE.match(line)
        if match is None:
            output.append(line)
            continue
        include_name = match.group(1)
        if Path(include_name).name == 'kernelParams.hpp':
            output.append('// kernelParams.hpp injected by mps_backend.py')
            continue
        include_path = _resolve_local_include(include_name, current_dir, search_dirs)
        if include_path in active_stack:
            chain = ' -> '.join(str(path) for path in (*active_stack, include_path))
            raise RuntimeError(f'Recursive Metal include detected: {chain}')
        output.append(f'// BEGIN INLINED INCLUDE: {include_path}')
        output.append(_inline_local_includes(
            include_path.read_text(encoding='utf-8'),
            current_dir=include_path.parent,
            search_dirs=search_dirs,
            active_stack=(*active_stack, include_path),
        ))
        output.append(f'// END INLINED INCLUDE: {include_path}')
    return '\n'.join(output)


def _find_kernel_params(root: Path, search_dirs: tuple[Path, ...]) -> Path:
    candidates = [
        root / 'kernelParams.hpp',
        root.parent / 'cpp/kernelParams.hpp',
        *(directory / 'kernelParams.hpp' for directory in search_dirs),
        Path(__file__).resolve().with_name('kernelParams.hpp'),
    ]
    for candidate in candidates:
        candidate = candidate.resolve()
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f'kernelParams.hpp was not found; searched: {candidates}')


def _assemble_metal_source(
    source_body: str,
    compiler_options: Iterable[Any],
    source_root: os.PathLike[str] | str,
) -> str:
    root = Path(source_root).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f'Metal source directory does not exist: {root}')
    search_dirs = (root, root.parent, root / 'include')
    expanded = _inline_local_includes(source_body, current_dir=root, search_dirs=search_dirs)
    kernel_params = _find_kernel_params(root, search_dirs).read_text(encoding='utf-8')
    if 'OMEGA_KERNEL_PARAMS_HPP_INCLUDED' not in kernel_params:
        kernel_params = (
            '#ifndef OMEGA_KERNEL_PARAMS_HPP_INCLUDED\n'
            '#define OMEGA_KERNEL_PARAMS_HPP_INCLUDED 1\n\n'
            + kernel_params
            + '\n#endif // OMEGA_KERNEL_PARAMS_HPP_INCLUDED\n'
        )
    return (
        _macro_preamble(compiler_options)
        + '#include <metal_stdlib>\nusing namespace metal;\n\n'
        + kernel_params
        + '\n'
        + expanded
        + '\n'
    )


def _compile_shader_cached(torch: Any, source: str, label: str) -> Any:
    digest = hashlib.sha256(source.encode('utf-8')).hexdigest()
    library = _SHADER_CACHE.get(digest)
    if library is None:
        try:
            library = torch.mps.compile_shader(source)
        except Exception as exc:
            raise RuntimeError(
                f'Metal {label} projector compilation failed. Source SHA-256: {digest}'
            ) from exc
        _SHADER_CACHE[digest] = library
    return library


def _without_compile_define(options: Iterable[Any], name: str) -> tuple[str, ...]:
    filtered: list[str] = []
    wanted = f'-D{name}'
    for option in options:
        tokens = [token for token in _iter_option_tokens((option,))
                  if token != wanted and not token.startswith(wanted + '=')]
        if tokens:
            filtered.append(' '.join(tokens))
    return tuple(filtered)


def _upload_static_buffers(self: Any, torch: Any) -> None:
    self.mps_empty_float32 = torch.empty(0, dtype=torch.float32, device='mps')
    self.mps_empty_uint8 = torch.empty(0, dtype=torch.uint8, device='mps')
    self.mps_empty_uint16 = torch.empty(0, dtype=torch.uint16, device='mps')
    self.mps_empty_uint32 = torch.empty(0, dtype=torch.uint32, device='mps')
    self.dummy_buffer = self.mps_empty_float32
    self.d_Sens = torch.zeros(1, dtype=torch.float32, device='mps')

    from .init import _initialize_coordinate_buffers, _initialize_detector_vector_buffers
    _initialize_coordinate_buffers(self, lambda value: _mps_tensor_from_numpy(torch, value, np.float32))
    _initialize_detector_vector_buffers(
        self,
        lambda value: _mps_tensor_from_numpy(torch, value, np.uint32),
        self.mps_empty_uint32,
    )
    if getattr(self, 'listmode', 0) > 0 and not getattr(self, 'useIndexBasedReconstruction', False) and not getattr(self, 'loadTOF', False):
        x = np.asarray(self.x, dtype=np.float32).ravel(order='F')
        self.d_x = [[self.mps_empty_float32] * self.subsets for _ in range(self.Nt)]
        for timestep in range(self.Nt):
            for subset in range(self.subsets):
                index = timestep * self.subsets + subset
                start = self.nMeas[index] * 6
                stop = self.nMeas[index + 1] * 6
                values = x[start:stop]
                if values.size:
                    self.d_x[timestep][subset] = _mps_tensor_from_numpy(torch, values, np.float32)

    self.d_rayShiftsDetector = _mps_tensor_from_numpy(torch, getattr(self, 'rayShiftsDetector', np.empty(0)), np.float32) if not _empty_value(getattr(self, 'rayShiftsDetector', np.empty(0))) else self.mps_empty_float32
    self.d_rayShiftsSource = _mps_tensor_from_numpy(torch, getattr(self, 'rayShiftsSource', np.empty(0)), np.float32) if not _empty_value(getattr(self, 'rayShiftsSource', np.empty(0))) else self.mps_empty_float32
    self.d_TOFCenter = _mps_tensor_from_numpy(torch, getattr(self, 'TOFCenter', np.empty(0)), np.float32) if not _empty_value(getattr(self, 'TOFCenter', np.empty(0))) else self.mps_empty_float32
    self.d_V = _mps_tensor_from_numpy(torch, getattr(self, 'V', np.empty(0)), np.float32) if not _empty_value(getattr(self, 'V', np.empty(0))) else self.mps_empty_float32

    attenuation = np.asarray(getattr(self, 'vaimennus', np.empty(0)), dtype=np.float32).ravel(order='F')
    measurement_attenuation = attenuation if self.attenuation_correction and not self.CTAttenuation else np.empty(0, dtype=np.float32)
    normalization = np.asarray(getattr(self, 'normalization', np.empty(0)), dtype=np.float32).ravel(order='F')
    corr_vector = np.asarray(getattr(self, 'corrVector', np.empty(0)), dtype=np.float32).ravel(order='F')
    offset_limit = np.asarray(getattr(self, 'OffsetLimit', np.empty(0)), dtype=np.float32).ravel(order='F')
    xy_index = np.asarray(getattr(self, 'xy_index', np.empty(0)), dtype=np.uint32).ravel(order='F')
    z_index = np.asarray(getattr(self, 'z_index', np.empty(0)), dtype=np.uint16).ravel(order='F')
    tr_index = np.asarray(getattr(self, 'trIndex', np.empty(0)), dtype=np.uint16).ravel(order='F')
    ax_index = np.asarray(getattr(self, 'axIndex', np.empty(0)), dtype=np.uint16).ravel(order='F')
    tof_index = np.asarray(getattr(self, 'TOFIndices', np.empty(0)), dtype=np.uint8).ravel(order='F')
    raw_data_length = np.asarray(getattr(self, 'LL', np.empty(0)), dtype=np.uint16).ravel(order='F')

    self.d_attenuation = [[self.mps_empty_float32] * self.subsets for _ in range(self.Nt)]
    self.d_norm = [[self.mps_empty_float32] * self.subsets for _ in range(self.Nt)]
    self.d_scatter = [[self.mps_empty_float32] * self.subsets for _ in range(self.Nt)]
    self.d_T = [[self.mps_empty_float32] * self.subsets for _ in range(self.Nt)]
    self.d_xyindex = [[self.mps_empty_uint32] * self.subsets for _ in range(self.Nt)]
    self.d_zindex = [[self.mps_empty_uint16] * self.subsets for _ in range(self.Nt)]
    self.d_trIndex = [[self.mps_empty_uint16] * self.subsets for _ in range(self.Nt)]
    self.d_axIndex = [[self.mps_empty_uint16] * self.subsets for _ in range(self.Nt)]
    self.d_TOFIndex = [[self.mps_empty_uint8] * self.subsets for _ in range(self.Nt)]
    self.d_L = [[self.mps_empty_uint16] * self.subsets for _ in range(self.Nt)]
    for timestep in range(self.Nt):
        for subset in range(self.subsets):
            index = timestep * self.subsets + subset
            measurement_start = self.nTotMeas[index]
            measurement_stop = self.nTotMeas[index + 1]
            projection_start = self.nMeas[index]
            projection_stop = self.nMeas[index + 1]
            if measurement_attenuation.size:
                self.d_attenuation[timestep][subset] = _mps_tensor_from_numpy(torch, measurement_attenuation[measurement_start:measurement_stop], np.float32)
            if normalization.size and getattr(self, 'normalization_correction', False):
                if getattr(self, 'SPECT', False) and int(getattr(self, 'normZ', 1)) == int(getattr(self, 'nHeads', 1)):
                    self.d_norm[timestep][subset] = _mps_tensor_from_numpy(torch, normalization, np.float32)
                else:
                    self.d_norm[timestep][subset] = _mps_tensor_from_numpy(torch, normalization[measurement_start:measurement_stop], np.float32)
            if corr_vector.size:
                self.d_scatter[timestep][subset] = _mps_tensor_from_numpy(torch, corr_vector[measurement_start:measurement_stop], np.float32)
            if offset_limit.size:
                self.d_T[timestep][subset] = _mps_tensor_from_numpy(torch, offset_limit[projection_start:projection_stop], np.float32)
            if xy_index.size:
                self.d_xyindex[timestep][subset] = _mps_tensor_from_numpy(torch, xy_index[projection_start:projection_stop], np.uint32)
            if z_index.size:
                self.d_zindex[timestep][subset] = _mps_tensor_from_numpy(torch, z_index[projection_start:projection_stop], np.uint16)
            if tr_index.size:
                self.d_trIndex[timestep][subset] = _mps_tensor_from_numpy(torch, tr_index[projection_start * 2:projection_stop * 2], np.uint16)
            if ax_index.size:
                self.d_axIndex[timestep][subset] = _mps_tensor_from_numpy(torch, ax_index[projection_start * 2:projection_stop * 2], np.uint16)
            if tof_index.size:
                self.d_TOFIndex[timestep][subset] = _mps_tensor_from_numpy(torch, tof_index[projection_start:projection_stop], np.uint8)
            if raw_data_length.size:
                self.d_L[timestep][subset] = _mps_tensor_from_numpy(torch, raw_data_length[projection_start:projection_stop], np.uint16)

    self.d_attenuation_image = self.mps_empty_float32
    if self.attenuation_correction and self.CTAttenuation and attenuation.size:
        self.d_attenuation_image = _mps_tensor_from_numpy(torch, attenuation, np.float32)

    mask_fp = getattr(self, 'maskFP', np.empty(0))
    self.d_maskFP = [[self.mps_empty_uint8] * self.subsets for _ in range(self.Nt)]
    if getattr(self, 'useMaskFP', False) and not _empty_value(mask_fp):
        mask_fp = np.asarray(mask_fp, dtype=np.uint8).ravel(order='F')
        frame_stride = int(getattr(self, 'measurement_nRowsD', self.nRowsD) * getattr(self, 'measurement_nColsD', self.nColsD))
        for timestep in range(self.Nt):
            for subset in range(self.subsets):
                if getattr(self, 'SPECT', False) and int(getattr(self, 'maskFPZ', 1)) == int(getattr(self, 'nHeads', 1)):
                    values = mask_fp
                elif int(getattr(self, 'maskFPZ', 1)) > 1:
                    index = timestep * self.subsets + subset
                    start = self.nMeas[index] * frame_stride
                    stop = self.nMeas[index + 1] * frame_stride
                    values = mask_fp[start:stop]
                else:
                    values = mask_fp
                self.d_maskFP[timestep][subset] = _mps_tensor_from_numpy(torch, values, np.uint8)
    self.d_maskBP = self.mps_empty_uint8
    if self.useMaskBP and self.maskBP.size:
        self.d_maskBP = _mps_tensor_from_numpy(torch, self.maskBP.ravel(order='F'), np.uint8)

    if int(getattr(self, 'FPType', 0)) == 6 or int(getattr(self, 'BPType', 0)) == 6:
        if isinstance(self.gFilter, (list, tuple)) and len(self.gFilter):
            self.d_gFilter = [
                _mps_tensor_from_numpy(torch, value, np.float32) for value in self.gFilter
            ]
        else:
            self.d_gFilter = _mps_tensor_from_numpy(torch, self.gFilter, np.float32)


def _geometry_buffer(self: Any, name: str, timestep: int, subset: int) -> Any:
    table = self.d_x if name == 'x' else self.d_z
    selected_subset = subset
    if name == 'x':
        subset_geometry = (getattr(self, 'CT', False) or getattr(self, 'SPECT', False)) and getattr(self, 'listmode', 0) == 0
        subset_geometry = subset_geometry or (getattr(self, 'listmode', 0) > 0 and not getattr(self, 'useIndexBasedReconstruction', False) and getattr(self, 'loadTOF', False))
        if not subset_geometry:
            selected_subset = 0
    else:
        subset_geometry = ((getattr(self, 'CT', False) or getattr(self, 'SPECT', False) or getattr(self, 'PET', False)) and getattr(self, 'listmode', 0) == 0)
        if not subset_geometry and not (getattr(self, 'listmode', 0) > 0 and not getattr(self, 'useIndexBasedReconstruction', False)):
            selected_subset = 0
    value = table[timestep][selected_subset]
    return self.mps_empty_float32 if value is None else value


def init_mps_projector(
    self: Any,
    *,
    source_root: os.PathLike[str] | str | None = None,
    source_fp: str | None = None,
    source_bp: str | None = None,
    options_fp: Iterable[Any] = (),
    options_bp: Iterable[Any] = (),
) -> None:
    import torch

    _validate_configuration(self)
    if not torch.backends.mps.is_available():
        raise RuntimeError('PyTorch MPS is not available on this machine.')
    if (int(self.FPType) != 6 or int(self.BPType) != 6) and not hasattr(torch.mps, 'compile_shader'):
        raise RuntimeError('This backend requires torch.mps.compile_shader().')

    self.useImages = False # PyTorch binds arrays as Metal buffers

    self.no_norm = 1
    self.mSize = int(getattr(self, 'measurement_nRowsD', self.nRowsD) * getattr(self, 'measurement_nColsD', self.nColsD) * self.nProjections)

    _upload_static_buffers(self, torch)
    # Type 6 uses the native MPS implementation for that direction.  Hybrid projectors still need a compiled shader for the other direction (for example type 61 is rotation FP + Siddon BP), so compile each side independently instead of returning early when FPType is 6.
    if (self.FPType != 6 or self.BPType != 6) and source_root is None:
        raise ValueError('Metal projector source is required for projector types 1-4')
    options_fp = _without_compile_define(tuple(options_fp), 'USEIMAGES')
    options_bp = _without_compile_define(tuple(options_bp), 'USEIMAGES')
    if self.FPType != 6:
        if source_fp is None:
            raise ValueError('Metal forward projector source is required for this configuration')
        complete_fp = _assemble_metal_source(source_fp, options_fp, source_root)
        self.mps_lib_fp = _compile_shader_cached(torch, complete_fp, 'forward')
        try:
            self.knlF = self.mps_lib_fp.projectorType123
        except AttributeError as exc:
            raise RuntimeError("Compiled Metal library does not expose 'projectorType123'.") from exc
    if self.BPType != 6:
        if source_bp is None:
            raise ValueError('Metal backward projector source is required for this configuration')
        complete_bp = _assemble_metal_source(source_bp, options_bp, source_root)
        self.mps_lib_bp = _compile_shader_cached(torch, complete_bp, 'backward')
        bp_name = 'projectorType123' if self.BPType in (1, 2, 3) else 'projectorType4Backward'
        try:
            self.knlB = getattr(self.mps_lib_bp, bp_name)
        except AttributeError as exc:
            raise RuntimeError(f"Compiled Metal library does not expose '{bp_name}'.") from exc

    if self.FPType != 6 or self.BPType != 6:
        self.d_scalar_params = []
        for timestep in range(self.Nt):
            per_subset = []
            for subset in range(self.subsets):
                per_volume = []
                for volume in range(int(self.nMultiVolumes) + 1):
                    packed = np.frombuffer(_pack_scalar_kernel_params(self, timestep, subset, volume), dtype=np.uint8).copy()
                    per_volume.append(torch.as_tensor(packed, device='mps'))
                per_subset.append(per_volume)
            self.d_scalar_params.append(per_subset)

def _kernel_args(
    self: Any,
    scalar_params: Any,
    dynamic_input: Any,
    output: Any,
    subset: int,
    timestep: int,
    direction: str,
) -> list[Any]:
    """Bind every Metal resource slot, using typed empty buffers when inactive."""
    empty_f = self.mps_empty_float32
    atten = self.d_attenuation_image if getattr(self, 'CTAttenuation', False) else self.d_attenuation[timestep][subset]
    if (direction == 'forward' and self.FPType in (1, 2, 3)) or (direction == 'backward' and self.BPType in (1, 2, 3)):
        args = [empty_f] * 22
        args[0] = scalar_params
        args[1] = self.d_rayShiftsDetector
        args[2] = self.d_rayShiftsSource
        args[3] = self.d_TOFCenter
        args[4] = self.d_V
        args[5] = atten
        args[6] = self.d_maskFP[timestep][subset]
        args[7] = self.d_maskBP
        args[8] = _geometry_buffer(self, 'x', timestep, subset)
        args[9] = _geometry_buffer(self, 'z', timestep, subset)
        args[10] = self.d_norm[timestep][subset]
        args[11] = self.d_scatter[timestep][subset]
        args[12] = self.d_Sens
        args[13] = self.d_xyindex[timestep][subset]
        args[14] = self.d_zindex[timestep][subset]
        args[15] = self.d_trIndex[timestep][subset]
        args[16] = self.d_axIndex[timestep][subset]
        args[17] = self.d_TOFIndex[timestep][subset]
        args[18] = self.d_L[timestep][subset]
        args[19] = dynamic_input
        args[20] = output
        args[21] = self.d_detectorVector[timestep][subset]
    elif direction == 'backward' and self.BPType == 4:
        args = [empty_f] * 10
        args[0] = scalar_params
        args[1] = self.d_T[timestep][subset]
        args[2] = dynamic_input
        args[3] = getattr(self, 'mps_fdk_angle', empty_f)
        args[4] = output
        args[5] = _geometry_buffer(self, 'x', timestep, subset)
        args[6] = _geometry_buffer(self, 'z', timestep, subset)
        args[7] = self.d_Sens
        args[8] = self.d_norm[timestep][subset]
        args[9] = self.d_maskBP
    return args


def _require_mps_float32_contiguous(tensor: Any, name: str) -> Any:
    import torch
    if not isinstance(tensor, torch.Tensor):
        raise TypeError(f'{name} must be a PyTorch tensor')
    if tensor.device.type != 'mps':
        raise ValueError(f"{name} must be on device='mps', got {tensor.device}")
    if tensor.dtype != torch.float32:
        raise TypeError(f'{name} must use torch.float32, got {tensor.dtype}')
    return tensor.contiguous()


def _projection_size(self: Any, timestep: int, subset: int) -> int:
    if self.subsetType > 7 or self.subsets == 1:
        return int(getattr(self, 'measurement_nRowsD', self.nRowsD) * getattr(self, 'measurement_nColsD', self.nColsD) * self.nProjSubset[timestep, subset])
    return int(self.nMeasSubset[timestep, subset])


def forward_projection_mps(self: Any, f: Any, subset: int, timestep: int) -> Any:
    import torch
    if not 0 <= timestep < self.Nt or not 0 <= subset < self.subsets:
        raise IndexError('timestep and subset must identify an initialized frame/subset pair')
    _validate_configuration(self)
    volume_count = self.nMultiVolumes + 1
    inputs = list(f) if isinstance(f, (list, tuple)) else [f] * volume_count
    if len(inputs) != volume_count:
        raise ValueError(f'Expected {volume_count} volume tensors, got {len(inputs)}')
    for volume in range(volume_count):
        image = _require_mps_float32_contiguous(inputs[volume], f'volume {volume}')
        if image.numel() != self.N[volume]:
            raise ValueError(f'Volume {volume} has {image.numel()} elements, expected {self.N[volume]}')
        inputs[volume] = image
    output = torch.zeros(_projection_size(self, timestep, subset), dtype=torch.float32, device='mps')
    for volume, image in enumerate(inputs):
        partial = torch.zeros_like(output)
        if self.FPType == 6:
            from .projfunctions import _type6_torch_ops, type6_forward
            with torch.no_grad(): # Disabling GradMode prevents MPS from retaining per-view state across OSEM subsets.
                type6_forward(self, image, partial, volume, subset, timestep, ops=_type6_torch_ops(self))
            torch.mps.synchronize()
        else:
            args = _kernel_args(self, self.d_scalar_params[timestep][subset][volume], image, partial, subset, timestep, 'forward')
            self.knlF(
                *args,
                threads=tuple(int(value) for value in self.globalSizeFP[timestep][subset]),
                group_size=tuple(int(value) for value in self.localSizeFP),
            )
        output += partial
    return output


def backward_projection_mps(self: Any, y: Any, subset: int, timestep: int) -> Any:
    import torch
    if not 0 <= timestep < self.Nt or not 0 <= subset < self.subsets:
        raise IndexError('timestep and subset must identify an initialized frame/subset pair')
    _validate_configuration(self)
    y = _require_mps_float32_contiguous(y, 'backprojection input')
    expected = _projection_size(self, timestep, subset)
    if y.numel() != expected:
        raise ValueError(f'Backprojection input has {y.numel()} elements, expected {expected}')
    outputs: list[Any] = []
    for volume in range(int(self.nMultiVolumes) + 1):
        output = torch.zeros(int(np.asarray(self.N).reshape(-1)[volume]), dtype=torch.float32, device='mps')
        if self.BPType == 6:
            from .projfunctions import _type6_torch_ops, type6_backward
            with torch.no_grad():
                type6_backward(self, y, output, volume, subset, timestep, ops=_type6_torch_ops(self))
            torch.mps.synchronize()
        else:
            args = _kernel_args(self, self.d_scalar_params[timestep][subset][volume], y, output, subset, timestep, 'backward')
            self.knlB(
                *args,
                threads=tuple(int(value) for value in self.globalSizeBP[timestep][subset][volume]),
                group_size=tuple(int(value) for value in self.localSizeBP),
            )
        outputs.append(output)
    return outputs[0] if self.nMultiVolumes == 0 else outputs
