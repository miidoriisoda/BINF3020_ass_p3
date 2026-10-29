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

    for base1, base2 in zip(reference_sequence, distant_sequence):
        # JC correction ignores indels
        if base1 == '-' or base2 == '-':
            continue

        valid_sites += 1
        if base1 != base2:
            mismatches += 1

    if valid_sites == 0:
        return 0.0

    # Calculate the proportion of differing sites (p-distance)
    p = mismatches / valid_sites

    # Jukes-Cantor is mathematically undefined for p >= 0.75 (saturation)
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
    raise NotImplementedError()
