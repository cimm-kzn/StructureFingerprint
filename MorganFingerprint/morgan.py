# -*- coding: utf-8 -*-
#
#  Copyright 2020 Aleksandr Sizov <murkyrussian@gmail.com>
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of MorganFingerprint.
#
#  MorganFingerprint is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from __future__ import annotations
from CGRtools.algorithms.morgan import tuple_hash
from collections import defaultdict, deque
from math import log2
from numpy import zeros
from sklearn.base import BaseEstimator, TransformerMixin
from typing import Any, Collection, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from CGRtools import MoleculeContainer


class MorganFingerprint(TransformerMixin, BaseEstimator):
    def __init__(self, radius: int = 4, length: int = 1024):
        """
        Utility tool for making Morgan-like fingerprints

        :param radius: maximum length of fragments
        :param length: bit string's length
        """
        self._radius = radius
        self._mask = length - 1
        self._length = length
        self._log = int(log2(length))

    def transform(self, x: Collection):
        fingerprints = zeros((len(x), self._length))

        for idx, mol in enumerate(x):
            fingerprint = fingerprints[idx]
            arr = self._fragments(self._bfs(mol), mol)
            hashes = {tuple_hash(tpl) for tpl in arr}

            for one in hashes:
                fingerprint[one & self._mask - 1] = 1
                fingerprint[one >> self._log & self._mask - 1] = 1

        return fingerprints

    def _bfs(self, molecule: MoleculeContainer) -> list[list[int]]:
        atoms = molecule._atoms
        bonds = molecule._bonds

        arr = [[x] for x in atoms]
        queue = deque(arr)
        while queue:
            now = queue.popleft()
            if len(now) >= self._radius:
                continue
            var = [now + [x] for x in bonds[now[-1]] if x not in now]
            arr.extend(var)
            queue.extend(var)
        return arr

    def _fragments(self, arr: list[list], molecule: MoleculeContainer) -> set[tuple[Union[int, Any], ...]]:
        atoms = {x: int(a) for x, a in molecule.atoms()}
        bonds = molecule._bonds
        cache = defaultdict(dict)
        out = set()
        for frag in arr:
            var = [atoms[frag[0]]]
            for x, y in zip(frag, frag[1:]):
                b = cache[x].get(y)
                if not b:
                    b = cache[x][y] = cache[y][x] = int(bonds[x][y])
                var.append(b)
                var.append(atoms[y])
            var = tuple(var)
            rev_var = var[::-1]
            if var > rev_var:
                out.add(var)
            else:
                out.add(rev_var)
        return out


__all__ = ['MorganFingerprint']
