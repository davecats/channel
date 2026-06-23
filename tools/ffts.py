from pathlib import Path
import sys

POSTPRO_DIR = Path(__file__).resolve().parents[1] / "postpro"
postpro_path = str(POSTPRO_DIR)
if postpro_path not in sys.path:
    sys.path.insert(0, postpro_path)

import numpy as np
import dask.array as da
import dask
import xarray as xr
from utils import mesh2k

from xrft.xrft import _new_dims_and_coords

def full_mesh2k(phi: float, node_size: int) -> np.ndarray:
    """Full Fourier wavenumber layout for padded inverse transforms."""
    if node_size % 2 == 0:
        half = node_size // 2
        k = np.empty(node_size)
        k[:half] = np.arange(0, half)*phi
        k[half:] = -np.arange(half, 0, -1)*phi
        return k
    return mesh2k(phi, (node_size - 1)//2)


"""Transform channel Fourier-space fields to physical space and back.

The code operates on dask-backed xr.DataArray objects and minimizes the number
of tasks while keeping operations chunk-local.
"""

def ifft_minimal(tmp: xr.DataArray, dim: str = "var") -> xr.DataArray:
    """
    Mimic xrft.ifft(tmp, dim=dim, true_amplitude=False, lag=0, shift=False)
    for a single dimension:
      - inverse FFT via np.fft.ifft
      - frequency coordinate via np.fft.fftfreq, as xrft._ifreq does
      - new dimension name freq_<dim>
    Assumes regular, increasing coord along `dim`.
    """
    axis = tmp.get_axis_num(dim)


    # Inverse FFT (no shift, no scaling)
    if tmp.data.chunks:
        f = dask.array.fft.ifft(tmp.data, axis=axis)
    else:
        f = np.fft.ifft(tmp.data, axis=axis)

    # Spacing from coord (xrft uses _diff_coord then abs(diff[0]))
    coord = tmp[dim].values
    diff = np.diff(coord)
    if diff.size == 0:
        raise ValueError(f"Need at least 2 points along {dim} to infer spacing")
    delta_x = float(abs(diff[0]))
    N = coord.size

    # Frequency coordinate: this is exactly what _ifreq gives for real=None, shift=False
    freq_coord = np.fft.fftfreq(N, delta_x)

    # New dim name: freq_<dim>, like xrft with prefix="freq_"
    freq_dim = f"freq_{dim}"

    # Coords: drop old Fourier coord, add freq_dim
    coords = {k: v for k, v in tmp.coords.items() if k != dim}
    coords[freq_dim] = freq_coord

    # Dims: replace dim with freq_dim
    new_dims = [freq_dim if d == dim else d for d in tmp.dims]

    return xr.DataArray(
        f,
        dims=new_dims,
        coords=coords,
        attrs=tmp.attrs,
        name=tmp.name,
    )

def layout_fourier_axis(
    data,
    axis: int,
    nodeSize: int,
    mirror: bool,
    nx: int = 0,
):


    def _pad_zeros_block(block, nodeSize_local, nx, mirror):
        if mirror:
            Nk = block.shape[-1]
            M = 2 * Nk - 1
            full = np.empty(block.shape[:-1] + (M,), dtype=block.dtype)

            conj_block = np.conjugate(block)
            full[..., : Nk - 1] = conj_block[..., Nk - 1:0:-1]
            full[..., Nk - 1:] = block[..., 0:]

            # Now apply a logical "roll" of -nx by explicit indexing:
            nx_mod = nx % M
            if nx_mod == 0:
                out = full
            else:
                out = np.empty_like(full)
                out[..., : M - nx_mod] = full[..., nx_mod:]
                out[..., M - nx_mod:] = full[..., :nx_mod]
            block=out

        M = block.shape[-1]
        nModes = (M - 1) // 2
        if 2 * nModes + 1 != M:
            raise ValueError(f"pad_zeros: expected odd length, got {M}")
        pad_len = nodeSize_local - 1 - 2 * nModes
        if pad_len < 0:
            raise ValueError("pad_zeros: nodeSize too small")

        left = block[..., : nModes + 1]     # (..., nModes+1)
        right = block[..., nModes + 1 :]    # (..., nModes)
        zeros = np.zeros(block.shape[:-1] + (pad_len,), dtype=block.dtype)
        return np.concatenate([left, zeros, right], axis=-1)

    size_axis = data.shape[axis]
    if len(data.chunks[axis]) != 1 or data.chunks[axis][0] != size_axis:
        data = data.rechunk({axis: size_axis})

    data_m = da.moveaxis(data, axis, -1)
    chunks = list(data_m.chunks)
    new_chunks = list(chunks)
    new_chunks[-1] = (nodeSize,)

    out_m = da.map_blocks(
        _pad_zeros_block,
        data_m,
        dtype=data_m.dtype,
        chunks=tuple(new_chunks),
        nodeSize_local=nodeSize,
        nx=nx,
        mirror=mirror
    )

    return da.moveaxis(out_m, -1, axis)

def my_ifft_Xarray(
    fourierField: xr.DataArray,
    dim: str,
    phi: float,
    nodeSize: int,
    out_dim: str | None = None,
    mirror: bool = False,
    nx = 0,
) -> xr.DataArray:
    """
    Inverse FFT along `dim`, optionally mirroring a folded spectrum first.

    mirror=False:
        - lays out the compact Fourier spectrum and pads unresolved modes.

    mirror=True:
        - first mirror+conj the folded spectrum along `dim`
        - then pads unresolved modes in the middle.
    """

    axis = fourierField.get_axis_num(dim)
    data = fourierField.data

    # layout: pad (and optionally mirror) along Fourier axis
    laid_out = layout_fourier_axis(
        data,
        axis=axis,
        nodeSize=nodeSize,
        mirror=mirror,
        nx=nx,
    )

    tmp = xr.DataArray(
        laid_out,
        dims=[d if d != dim else "var" for d in fourierField.dims],
        coords={k: v for k, v in fourierField.coords.items() if k != dim},
        attrs=fourierField.attrs,
        name=fourierField.name,
    )

    # assign 'var' coord via mesh2k (identical to old)
    var_coord = full_mesh2k(phi, nodeSize)
    if var_coord.size != nodeSize:
        raise ValueError(
            f"mesh2k returned {var_coord.size}, expected {nodeSize}"
        )
    tmp = tmp.assign_coords(var=var_coord)

    # IFFT
    physicalField = ifft_minimal(tmp, dim="var")
    physicalField = physicalField * nodeSize

    # rename freq_var -> physical dim
    if out_dim is None:
        new_dim = dim[-1]
    else:
        new_dim = out_dim

    physicalField = physicalField.rename({"freq_var": new_dim})
    physicalField = physicalField.assign_coords(
        {new_dim: np.arange(nodeSize) / nodeSize * (2 * np.pi / phi)}
    )

    return physicalField

def physical(xArray: xr.DataArray, dns: dict) -> xr.DataArray:
    # 1) kz IFFT (no mirror), same as old
    tmp = my_ifft_Xarray(
        xArray,
        dim="kz",
        phi=dns["beta0"],
        nodeSize=dns["nzd"],
        out_dim="z",
        mirror=False,   # folded kz is already handled by the spectrum layout
    )

    # 3) kx IFFT (no mirror), same as old
    physicalArray = my_ifft_Xarray(
        tmp,
        dim="kx_folded",
        phi=dns["alfa0"],
        nodeSize=2*dns["nxd"],
        out_dim="x",
        mirror=True,
        nx = dns["nx"]
    )

    return physicalArray.real


def _ifftshift_block(block, axis=-1):
    return np.fft.ifftshift(block, axes=axis)

def ifftshift_chunkwise(d: xr.DataArray, dim: str) -> xr.DataArray:
    """
    Apply np.fft.ifftshift along `dim` *within each chunk* of da.data.

    Does not merge chunks along dim; each chunk is independently shifted.
    """
    axis = d.get_axis_num(dim)
    data = d.data

    if isinstance(data, da.Array):
        # move axis to last position so block function is simple
        data_m = da.moveaxis(data, axis, -1)

        # apply blockwise ifftshift along last axis
        out_m = da.map_blocks(
            _ifftshift_block,
            data_m,
            dtype=data_m.dtype,
            axis=-1,
        )

        out = da.moveaxis(out_m, -1, axis)
    else:
        # pure NumPy
        out = np.fft.ifftshift(data, axes=axis)

    return xr.DataArray(
        out,
        dims=d.dims,
        coords=d.coords,
        attrs=d.attrs,
        name=d.name,
    )

def fft_minimal(
    da,
    dim=None,
    prefix="freq_",
    true_phase=True,
):

    axis_num = da.get_axis_num(dim)

    # verify even spacing of input coordinates
    diff = np.diff(da[dim])
    delta_x = np.abs(diff[0])
    lag_x = da[dim].data[len(da[dim].data) // 2]
    if delta_x == 0.0:
        raise ValueError(
            "Can't take Fourier transform because spacing in coordinate %s is zero"
            % dim
        )
    if true_phase:
        shifted = ifftshift_chunkwise(da, dim)
        f = dask.array.fft.fftn(
            shifted,
            axes=[axis_num],
        )
    else:
        f = dask.array.fft.fftn(da, axes=[axis_num])

    k = [np.fft.fftfreq(da.shape[axis_num], delta_x)]

    newcoords, swap_dims = _new_dims_and_coords(da, [dim], k, prefix)
    daft = xr.DataArray(
        f, dims=da.dims, coords=dict([c for c in da.coords.items() if c[0] != dim])
    )
    daft = daft.swap_dims(swap_dims).assign_coords(newcoords)

    up_dim = daft.dims[da.get_axis_num(dim)]
    if true_phase:
        daft = daft * xr.DataArray(
            np.exp(-1j * 2.0 * np.pi * newcoords[up_dim] * lag_x),
            dims=up_dim,
            coords={up_dim: newcoords[up_dim]},
        )
    daft[up_dim].attrs.update({"direct_lag": lag_x})

    return daft

def _compact_modes_block(block, modeSize):
    """
    Given a block with last axis length N, return a compact spectrum of length
    2*modeSize+1:

        out = [ block[..., 0:modeSize+1], block[..., N-modeSize:N] ]
    """
    N = block.shape[-1]
    if modeSize >= N:
        raise ValueError(f"modeSize {modeSize} must be < N {N}")
    left = block[..., :modeSize+1]     # shape (..., modeSize+1)
    right = block[..., N-modeSize:]    # shape (..., modeSize)
    return np.concatenate([left, right], axis=-1)

def compact_fourier_axis(
    daft: xr.DataArray,
    dim: str,
    modeSize: int,
    out_dim: str = "var",
) -> xr.DataArray:
    """
    Take xrft-style fft output along `dim` (length N),
    and build a compact symmetric spectrum of length 2*modeSize+1
    by taking [0..modeSize] and the last `modeSize` entries.

    Returns a new DataArray with dimension name `out_dim` along that axis.
    """
    axis = daft.get_axis_num(dim)
    data = daft.data

    # move transform axis to last
    if isinstance(data, da.Array):
        data_m = da.moveaxis(data, axis, -1)
        N = data_m.shape[-1]
        # ensure single chunk along that axis (to keep the kernel simple)
        if len(data_m.chunks[-1]) != 1 or data_m.chunks[-1][0] != N:
            data_m = data_m.rechunk({data_m.ndim - 1: N})

        out_len = 2 * modeSize + 1
        out_chunks = list(data_m.chunks)
        out_chunks[-1] = (out_len,)

        out_m = da.map_blocks(
            _compact_modes_block,
            data_m,
            dtype=data_m.dtype,
            chunks=tuple(out_chunks),
            modeSize=modeSize,
        )
        out_data = da.moveaxis(out_m, -1, axis)
    else:
        arr = np.moveaxis(np.asarray(data), axis, -1)
        out_arr = _compact_modes_block(arr, modeSize)
        out_data = np.moveaxis(out_arr, -1, axis)

    # build new dims/coords
    new_dims = [out_dim if d == dim else d for d in daft.dims]
    coords = {name: c for name, c in daft.coords.items() if name != dim}
    # we assign the logical k-grid later via mesh2k, so no coord for out_dim here

    return xr.DataArray(
        out_data,
        dims=new_dims,
        coords=coords,
        attrs=daft.attrs,
        name=daft.name,
    )

def my_fft_Xarray(
    physicalField: xr.DataArray,
    dim: str,
    phi: float,
    modeSize: int,
) -> xr.DataArray:
    nodeSize = physicalField.sizes[dim]

    # forward FFT (your optimized version without true_phase)
    fourierField = fft_minimal(physicalField, dim=dim, prefix="freq_", true_phase=True)
    fourierField = fourierField / nodeSize

    freq_dim = f"freq_{dim}"
    # compact spectrum along freq_dim into length 2*modeSize+1
    compact = compact_fourier_axis(
        fourierField,
        dim=freq_dim,
        modeSize=modeSize,
        out_dim="var",
    )

    # final naming and coords as before
    out = compact.rename({"var": "k" + dim})
    out = out.assign_coords({"k" + dim: mesh2k(phi, modeSize)})

    return out

def fourier(xArray, dns):
    tmpField1 = my_fft_Xarray(xArray, 'x', dns['alfa0'], dns['nx'])
    tmpField2 = tmpField1.isel(kx = slice(None,dns['nx']+1)).rename({'kx': 'kx_folded'})

    fourierArray = my_fft_Xarray(tmpField2, 'z', dns['beta0'], dns['nz'])

    return fourierArray
