def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    import math
    # First check sequences are the same length
    if len(reference_sequence) != len(distant_sequence):
        raise ValueError("Sequences must be aligned and of equal length.")

    # Count the number of mismatches
    mismatches = 0
    valid_sites = 0
    VALID = {"A", "G", "C", "T"}

    ref_seq = reference_sequence.upper()
    dist_seq = distant_sequence.upper()
    for base1, base2 in zip(ref_seq, dist_seq):
        # JC correction ignores indels
        if base1 == '-' or base2 == '-':
            continue
        elif base1 not in VALID or base2 not in VALID:
            continue
        valid_sites += 1
        if base1 != base2:
            mismatches += 1

    if valid_sites == 0:
        return 0.0

    # Calculate proportion of differing sites
    p = mismatches / valid_sites

    # JC undefined for p >= 0.75
    if p >= 0.75:
        return float('inf')

    # Apply the Jukes-Cantor correction formula
    distance = -0.75 * math.log(1 - (4 / 3) * p)
    return distance * valid_sites


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    import math
    VALID = {"A", "G", "C", "T"}
    PURINES = {"A", "G"}
    transitions = 0
    transversions = 0
    valid = 0

    ref_seq = reference_sequence.upper()
    dist_seq = distant_sequence.upper()
    for a, b in zip(ref_seq, dist_seq):
        if a not in VALID or b not in VALID:
            continue
        valid += 1
        if a == b:
            continue
        if (a in PURINES) == (b in PURINES):
            # purine to purine / pyrimidine to pyrimidine
            transitions += 1
        else:
            transversions += 1
    if valid == 0:
        return 0.0
    P = transitions / valid
    Q = transversions / valid
    arg1 = 1 - 2 * P - Q
    arg2 = 1 - 2 * Q
    if arg1 <= 0 or arg2 <= 0:
        return float("inf")
    d = -0.5 * math.log(arg1) - 0.25 * math.log(arg2)
    return d * valid
