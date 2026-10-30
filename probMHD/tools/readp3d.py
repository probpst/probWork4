#!/usr/bin/env python3
"""Readers for the probMHD PLOT3D files (Fortran unformatted, little-endian, float64)
and the monitor / slice files.

    from readp3d import read_grid, read_function, read_names, read_monitor
"""
import numpy as np


def _records(path):
    buf = open(path, 'rb').read()
    pos = 0
    while pos < len(buf):
        n = int(np.frombuffer(buf, np.int32, 1, pos)[0]); pos += 4
        yield buf[pos:pos+n]; pos += n
        n2 = int(np.frombuffer(buf, np.int32, 1, pos)[0]); pos += 4
        if n != n2:
            raise ValueError('record markers disagree')


def read_grid(path):
    rec = _records(path)
    ni, nj, nk = np.frombuffer(next(rec), np.int32, 3)
    raw = np.frombuffer(next(rec), np.float64)
    n = ni*nj*nk
    X, Y, Z = (raw[i*n:(i+1)*n].reshape((ni, nj, nk), order='F') for i in range(3))
    return X, Y, Z


def read_function(path):
    rec = _records(path)
    ni, nj, nk, nv = np.frombuffer(next(rec), np.int32, 4)
    raw = np.frombuffer(next(rec), np.float64)
    n = ni*nj*nk
    return [raw[i*n:(i+1)*n].reshape((ni, nj, nk), order='F') for i in range(nv)]


def read_names(path):
    return [l.strip() for l in open(path) if l.strip()]


def read_monitor(path):
    d = np.loadtxt(path)
    cols = ['step', 't', 'dt', 'cfl', 'ekin', 'emag', 'eth', 'etot', 'mass', 'vmax', 'divb', 'emagpert', 'tA', 'rhomin', 'pmin']
    return {c: d[:, i] for i, c in enumerate(cols) if i < d.shape[1]}


def read_slice(path):
    d = np.loadtxt(path)
    cols = ['i', 'j', 'rho', 'u', 'v', 'w', 'p', 'T', 'Bx', 'By', 'Bz', 'Jz', 'divB']
    return {c: d[:, k] for k, c in enumerate(cols)}
