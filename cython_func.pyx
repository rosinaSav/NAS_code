#compile using python3 cython_setup.py build_ext --inplace

import numpy as np

def calc_density_c(motif, sequence, double motif_length, return_count = False):
    cdef float density, match_length, seq_length
    if not sequence:
        return(None, None)
    matches = np.array([np.arange(hit.start(), hit.start() + motif_length) for hit in motif.finditer(sequence)])
    matches = np.unique(matches)
    match_length = len(matches)
    if return_count:
        return(match_length, matches)
    seq_length = len(sequence)
    density = match_length/seq_length
    return(density, matches)

def calc_density_for_concat_c(motif, sequence, double motif_length):
    cdef float density, match_length, seq_length
    if not sequence:
        return(None)
    matches = np.array([np.arange(hit.start(), hit.start() + motif_length) for hit in motif.finditer(sequence)])
    return(matches)

def calc_density_for_concat_several_c(motifs, sequence, motif_lengths):
    cdef float density, match_length, seq_length
    if not sequence:
        return(None)
    matches = [np.array([np.arange(hit.start(), hit.start() + motif_lengths[j]) for hit in motifs[j].finditer(sequence)]) for j in range(len(motifs))]
    return(matches)
