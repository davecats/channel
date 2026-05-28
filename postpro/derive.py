import numpy as np
import xarray as xr
import dask.array as da
from numba import njit
from scipy.linalg import solve_banded
from scipy.linalg import solve as sp_solve

"""Banded vector-matrix multiplication and derivatives"""

DERIVS = np.dtype([('d0', np.float64), ('d1', np.float64), ('d2', np.float64), ('d4', np.float64)]) 
def setup_derivatives(dns,y):                                                              
    ny = dns['ny']                                                              
    der = np.zeros([ny+3,ny+3],DERIVS)
    matder=np.zeros([5,5],np.float64)
    rhsder=np.zeros([5],np.float64)  

    for iy in range(2,ny+1):
        matder = np.fromfunction(lambda i,j: (y[iy-2+j]-y[iy])**(4.0-np.float64(i)), (5,5), dtype=int)
        rhsder *= 0
        rhsder[0] = np.float64(24)
        der[iy,iy-2:iy+3]['d4'] = sp_solve(matder,rhsder)
        matder = np.fromfunction(lambda i,j: (5.0-i)*(6.0-i)*(7.0-i)*(8.0-i)*((y[iy-2+j]-y[iy])**(4.0-np.float64(i))), (5,5), dtype=int)
        rhsder *= 0

        for i in range(5):
            rhsder[i] = np.sum(der[iy,iy-2:iy+3]['d4']*(y[iy-2:iy+3]-y[iy])**np.float64(8-i))
            
        der[iy,iy-2:iy+3]['d0'] = sp_solve(matder,rhsder)
        matder = np.fromfunction(lambda i,j: (y[iy-2+j]-y[iy])**(4-np.float64(i)), (5,5), dtype=int)
        rhsder *= 0

        for i in range(3):
            rhsder[i] = np.sum((4.0-i)*(3.0-i)*der[iy,iy-2:iy+3]['d0']*(y[iy-2:iy+3]-y[iy])**np.float64(2-i))

        der[iy,iy-2:iy+3]['d2'] = sp_solve(matder,rhsder)
        rhsder *= 0

        for i in range(4):
            rhsder[i] = np.sum((4.0-i)*der[iy,iy-2:iy+3]['d0']*(y[iy-2:iy+3]-y[iy])**np.float64(3-i))

        der[iy,iy-2:iy+3]['d1'] = sp_solve(matder,rhsder)
        
    # Bottom wall
    der[1,1]['d0'] = np.float64(1.0)
    der[0,0]['d0'] = np.float64(1.0)

    matder = np.fromfunction(lambda i,j: (y[0+j]-y[1])**(4-np.float64(i)), (5,5), dtype=int)
    rhsder = np.zeros([5],np.float64)
    rhsder[3] = np.float64(1.0)
    der[1,0:5]['d1'] = sp_solve(matder,rhsder)
    rhsder = np.zeros([5],np.float64)
    rhsder[2] = np.float64(2.0)
    der[1,0:5]['d2'] = sp_solve(matder,rhsder)

    matder = np.fromfunction(lambda i,j: (y[0+j]-y[0])**(4-np.float64(i)), (5,5), dtype=int)
    rhsder = np.zeros([5],np.float64)
    rhsder[3] = np.float64(1.0)
    der[0,0:5]['d1'] = sp_solve(matder,rhsder)
    rhsder = np.zeros([5],np.float64)
    rhsder[2] = np.float64(2.0)
    der[0,0:5]['d2'] = sp_solve(matder,rhsder)

    # Top wall
    der[ny+1,ny+1]['d0'] = np.float64(1.0)
    der[ny+2,ny+2]['d0'] = np.float64(1.0)

    matder = np.fromfunction(lambda i,j: (y[ny-2+j]-y[ny+1])**(4-np.float64(i)), (5,5), dtype=int)
    rhsder = np.zeros([5],np.float64)
    rhsder[3] = np.float64(1.0)
    der[ny+1,ny-2:]['d1'] = sp_solve(matder,rhsder)
    rhsder = np.zeros([5],np.float64)
    rhsder[2] = np.float64(2.0)
    der[ny+1,ny-2:]['d2'] = sp_solve(matder,rhsder)

    matder = np.fromfunction(lambda i,j: (y[ny-2+j]-y[ny+2])**(4-np.float64(i)), (5,5), dtype=int)
    rhsder = np.zeros([5],np.float64)
    rhsder[3] = np.float64(1.0)
    der[ny+2,ny-2:]['d1'] = sp_solve(matder,rhsder)
    rhsder = np.zeros([5],np.float64)
    rhsder[2] = np.float64(2.0)
    der[ny+2,ny-2:]['d2'] = sp_solve(matder,rhsder)

    return der

def get_banded_derivative_operators(config, y):
    assert isinstance(y, xr.DataArray)
    size = [4 if key != "d0" else 2 for key in ["d0", "d1", "d2"]]

    der = setup_derivatives(config, y.data)
    return [xr.DataArray(
        build_D_banded_from_dense(der[key], l=s, u=s),
        dims=["band", "y"],
        coords={"band": np.arange(1+2*s), "y": y.data})
        for key, s in zip(["d0", "d1", "d2"], size)]

def build_D_banded_from_dense(D_dense, l=2, u=2):
    n = D_dense.shape[0]
    ab = np.zeros((l + u + 1, n), dtype=D_dense.dtype)
    for j in range(n):
        i_min = max(0, j - u)
        i_max = min(n, j + l + 1)
        for i in range(i_min, i_max):
            k = u + i - j
            ab[k, j] = D_dense[i, j]
    return ab


@njit
def _banded_mv_core_numba(ab, v, l, u):
    """
    Banded matrix-vector multiplication:
      ab: (l+u+1, ny)
      v : (ny, K)
    ab storage: ab[u + i - j, j] = A[i,j].
    Returns out: (ny, K) = A @ v.
    """
    ny, K = v.shape
    out = np.zeros((ny, K), dtype=v.dtype)
    for j in range(ny):
        i_min = max(0, j - u)
        i_max = min(ny, j + l + 1)
        for i in range(i_min, i_max):
            k = u + i - j
            aij = ab[k, j]
            if aij != 0.0:
                for c in range(K):
                    out[i, c] += aij * v[j, c]
    return out
def fourier_dy_fast(
    field: xr.DataArray,
    Dnb_xr: xr.DataArray,
    d0_banded_xr: xr.DataArray,
    band_dim: str = "band",
    y_dim: str = "y",
    l_Dnb: int = 4,
    u_Dnb: int = 4,   # Dnb band size = 9
    l_d0: int = 2,
    u_d0: int = 2,    # d0 band size = 5
) -> xr.DataArray:
    """
    Compute dy via:
        rhs = 
        d0_banded_xr @ out = Dnb_xr @ field
    in a single dask.map_blocks over `field.data`.

    Dnb_xr       : banded differentiation operator, dims (band_dim, y_dim),
                   size(band_dim) == l_Dnb + u_Dnb + 1, unchunked
    d0_banded_xr : banded solver operator, dims (band_dim, y_dim),
                   size(band_dim) == l_d0 + u_d0 + 1, unchunked
    field        : DataArray with dim y_dim and any batch dims
    """

    # --- validate operators (no alignment along band_dim) ---
    def _check_op(op: xr.DataArray, name: str, l: int, u: int):
        
        assert list(op.dims) == [band_dim, y_dim], f"{name} must have only dims {band_dim!r} and {y_dim!r}"
        assert op.sizes[band_dim] == l + u + 1, (
            f"{name}.{band_dim!r} size must be l+u+1={l+u+1}, "
             f"found {op.sizes[band_dim]}")
        # ensure order is (band_dim, y_dim)
        assert isinstance(op.data, np.ndarray), f"{name} must be a single unchunked block (numpy array), got {type(op.data)}"

    _check_op(Dnb_xr, "Dnb_xr", l_Dnb, u_Dnb)
    _check_op(d0_banded_xr, "d0_banded_xr", l_d0, u_d0)

    assert y_dim in field.dims, f"{y_dim!r} must be in field.dims"

    # Check that the y coordinates match
    xr.align(Dnb_xr.y, d0_banded_xr.y, field.y, join="exact")

    y_axis = field.get_axis_num(y_dim)

    assert isinstance(field.data, da.Array), "field.data must be a dask array for map_blocks"
    y_chunks = field.data.chunks[y_axis]
    assert (
        len(y_chunks) == 1 and y_chunks[0] == field.sizes[y_dim]
    ), f"{y_dim} must have a single full-length chunk for map_blocks; got chunks {y_chunks}"

    def _block(field_block, Dnb_ab, d0_ab, l_Dnb, u_Dnb, l_d0, u_d0):
        """
        field_block : numpy array, block of field with y at axis y_axis
        Dnb_ab      : numpy array, banded differentiation operator
        d0_ab       : numpy array, banded solver operator
        """
        # move y axis to front and flatten batch dims
        v = np.moveaxis(field_block, y_axis, 0)  # (ny, *batch)
        ny = v.shape[0]
        batch_shape = v.shape
        v2 = v.reshape(ny, -1)  # (ny, K)

        # 1) rhs = Dnb_ab @ field
        rhs2 = _banded_mv_core_numba(Dnb_ab, v2, l_Dnb, u_Dnb)

        # 2) out2 = d0^{-1} @ rhs via solve_banded
        out2 = solve_banded((l_d0, u_d0), d0_ab, rhs2)

        # restore shape and axes
        out = out2.reshape(batch_shape)
        out = np.moveaxis(out, 0, y_axis)
        return out

    out_data = da.map_blocks(
        _block,
        field.data,
        Dnb_xr.data,
        d0_banded_xr.data,
        l_Dnb,
        u_Dnb,
        l_d0,
        u_d0,
        dtype=field.data.dtype,
    )

    out = xr.DataArray(
        out_data,
        dims=field.dims,
        coords=field.coords,
        name=field.name,
        attrs=field.attrs,
    )
    return out