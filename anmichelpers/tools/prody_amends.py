"""
Created on 03/06/2015

Some canibalized functions from ProDy, slightly modified to work as I want :D

@author: vgil
"""
from prody.atomic.residue import Residue
from prody.measure.measure import getDihedral, calcDistance

def calcPhi(residue, radian=False, dist=4.1):
    if not isinstance(residue, Residue):
        raise TypeError('{0} must be a Residue instance')

    C_, N, CA, C = getPhiAtoms(residue, dist=dist)

    return getDihedral(C_._getCoords(), N._getCoords(), CA._getCoords(),
                       C._getCoords(), radian)


def getPhiAtoms(residue, dist=4.1):

    prev = residue.getPrev()
    try:
        isaa = prev.isaminoacid
    except AttributeError:
        raise ValueError('{0} is a terminal residue'.format(str(residue)))

    C_ = prev['C']
    if C_ is None:
        raise ValueError('{0} does not have C atom'.format(str(prev)))
    N = residue['N']
    if N is None:
        raise ValueError('{0} does not have N atom'.format(str(residue)))
    CA = residue['CA']
    if CA is None:
        raise ValueError('{0} does not have CA atom'.format(str(residue)))
    C = residue['C']
    if C is None:
        raise ValueError('{0} does not have C atom'.format(str(residue)))
    CA_ = prev['CA']
    if CA_ is None:
        raise ValueError('{0} does not have CA atom'.format(str(prev)))

    if dist and dist < calcDistance(CA, CA_):
        raise ValueError('{0} and {1} does not seem to be connected'
                         .format(str(residue), str(prev)))

    return C_, N, CA, C


def calcPsi(residue, radian=False, dist=4.1):

    if not isinstance(residue, Residue):
        raise TypeError('{0} must be a Residue instance')

    N, CA, C, _N = getPsiAtoms(residue, dist = dist)

    return getDihedral(N._getCoords(), CA._getCoords(), C._getCoords(),
                       _N._getCoords(), radian)


def getPsiAtoms(residue, dist=4.1):

    N = residue['N']
    if N is None:
        raise ValueError('{0} does not have N atom'.format(str(residue)))
    CA = residue['CA']
    if CA is None:
        raise ValueError('{0} does not have CA atom'.format(str(residue)))
    C = residue['C']
    if C is None:
        raise ValueError('{0} does not have C atom'.format(str(residue)))

    next_res = residue.getNext()
    if next_res is not None:
        _N = next_res['N']
        CA_ = next_res['CA']
        if CA_ is None:
            raise ValueError('{0} does not have CA atom'.format(str(next_res)))
        if dist and dist < calcDistance(CA, CA_):
            raise ValueError('{0} and {1} does not seem to be connected'
                         .format(str(residue), str(next_res)))
    else:    
        OXT = residue['OXT']
        if OXT is None:
            raise ValueError('{0} terminal residue has not OXT atom'.format(str(next_res)))
        else:
            _N = OXT

    return N, CA, C, _N
