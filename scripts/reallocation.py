"""Reproduce the upstream Fork-form operator; aggregate by native-cell area.

Upstream modules are read from an explicit checkout, without modifying it.
The large internal projected grid is a numerical workspace, not a fine-scale
outage product. Source-cell equivalence is tested, never enforced by rescaling.
"""
import importlib
import math
from pathlib import Path
import sys

import numpy as np
from pyproj import Transformer
from rasterio.transform import from_origin
from scipy import sparse
import shapely


def upstream_modules(path):
    src = Path(path).resolve() / 'src'
    if not (src / 'nocturne/disaggregate/operator.py').exists():
        raise FileNotFoundError('Expected the original ntl-psf-disaggregation checkout')
    # Do not create __pycache__ or any output in the original work repository.
    old = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(src))
    try:
        modules = {key: importlib.import_module(f'nocturne.disaggregate.{key}')
                   for key in ('operator', 'kernels', 'overture_bundle')}
        for module in modules.values():
            if not Path(module.__file__).resolve().is_relative_to(src):
                raise ValueError('A different nocturne installation was already imported')
    finally:
        sys.dont_write_bytecode = old
        sys.path.remove(str(src))
    return modules


def native_overlap_matrix(cells, fine_transform, fine_shape, target_crs='EPSG:32616', progress=False):
    """Exact projected polygon/pixel intersection weights for native-cell means.

    Geographic edges are densified before projection. Matrix row sums must be
    one: a partially covered parent cell is not silently renormalized.
    """
    transformer = Transformer.from_crs('EPSG:4326', target_crs, always_xy=True)
    height, width = fine_shape
    resolution = fine_transform.a
    weights, indices, pointers, areas = [], [], [0], []
    for i, cell in enumerate(cells.itertuples(index=False)):
        left = -90 + cell.source_col / 240
        top = 40 - cell.source_row / 240
        polygon = shapely.segmentize(shapely.box(left, top-1/240, left+1/240, top), .0001)
        polygon = shapely.transform(polygon, transformer.transform, interleaved=False)
        xmin, ymin, xmax, ymax = polygon.bounds
        c0 = max(0, math.floor((xmin-fine_transform.c)/resolution))
        c1 = min(width, math.ceil((xmax-fine_transform.c)/resolution))
        r0 = max(0, math.floor((fine_transform.f-ymax)/resolution))
        r1 = min(height, math.ceil((fine_transform.f-ymin)/resolution))
        rr, cc = np.meshgrid(np.arange(r0,r1), np.arange(c0,c1), indexing='ij')
        rr, cc = rr.ravel(), cc.ravel()
        x = fine_transform.c + cc * resolution
        y = fine_transform.f - rr * resolution
        boxes = shapely.box(x, y-resolution, x+resolution, y)
        inside = shapely.covers(polygon, boxes)
        overlap = np.full(len(boxes), resolution**2)
        overlap[~inside] = shapely.area(shapely.intersection(boxes[~inside], polygon))
        keep = overlap > 1e-8
        area = polygon.area
        if not np.isclose(overlap[keep].sum(), area, rtol=1e-8, atol=1e-5):
            raise ValueError(f'Fine grid does not cover native cell {i}')
        weights.append(overlap[keep] / area)
        indices.append((rr[keep]*width+cc[keep]).astype('int32'))
        pointers.append(pointers[-1]+int(keep.sum()))
        areas.append(area)
        if progress and (i+1)%1000 == 0:
            print(f'Native polygon overlaps: {i+1}/{len(cells)}', flush=True)
    matrix = sparse.csr_matrix((np.concatenate(weights), np.concatenate(indices), np.array(pointers)),
                               shape=(len(cells), height*width))
    if not np.allclose(np.asarray(matrix.sum(axis=1)).ravel(), 1, rtol=0, atol=1e-8):
        raise AssertionError('Native overlap weights do not sum to one')
    return matrix, np.array(areas)


def aggregate_complete(matrix, values, tolerance=1e-6):
    """Area-weighted native means, requiring effectively complete fine support."""
    v = np.asarray(values).ravel()
    support = np.asarray(matrix @ np.isfinite(v).astype('float64')).ravel()
    mean = np.asarray(matrix @ np.nan_to_num(v, nan=0.)).ravel()
    mean[support < 1-tolerance] = np.nan
    return mean, support


def fine_source_indices(transform, shape, source_row_start, source_col_start, source_shape):
    """Nearest-source-cell assignment at projected fine-pixel centers."""
    to_geo = Transformer.from_crs('EPSG:32616', 'EPSG:4326', always_xy=True)
    height, width = shape
    result = np.empty(shape, dtype='int32')
    for r0 in range(0,height,128):
        rr, cc = np.meshgrid(np.arange(r0,min(r0+128,height)), np.arange(width), indexing='ij')
        x = transform.c+(cc+.5)*transform.a
        y = transform.f+(rr+.5)*transform.e
        lon, lat = to_geo.transform(x,y)
        rows = np.floor((40-lat)*240).astype(int)-source_row_start
        cols = np.floor((lon+90)*240).astype(int)-source_col_start
        if np.any((rows<0)|(rows>=source_shape[0])|(cols<0)|(cols>=source_shape[1])):
            raise ValueError('Expanded raw VIIRS crop does not cover the working grid')
        result[r0:r0+len(rr)] = rows*source_shape[1]+cols
    return result


def allocate_nightly(source, gain, kernel, correlate, geometric_support, tolerance=1e-6):
    """Algebraically identical allocation component of upstream Fork-form.

    Static proxy denominator/gain and geometric support are precomputed once.
    Neighborhood support still depends on EACH night's strict valid radiance.
    """
    valid = np.isfinite(source)
    if not valid.any():
        empty = np.full(source.shape, np.nan, dtype='float32')
        return empty, empty.copy()
    support = correlate(valid.astype('float32'), kernel.weights)
    weighted_sum = correlate(np.where(valid,source,0).astype('float32'), kernel.weights)
    mean = np.full(source.shape,np.nan,dtype='float32')
    allowed = (support>=1-tolerance)&(geometric_support>=1-tolerance)&np.isfinite(gain)
    np.divide(weighted_sum,support,out=mean,where=allowed)
    allocation = np.maximum(gain * mean,0).astype('float32')
    return allocation, mean
