"""
Utilities to combine two genotype readers.
"""

from typing import Iterable, Iterator, Tuple, Callable
from multiprocessing import Process, Queue

from .core import Genotypes, Variant


ChromComparator = Callable[[str, str], int]


def chrom_to_int(chrom: str) -> int:
    try:
        return int(chrom)
    except ValueError:
        if chrom == "X":
            return 23
        if chrom == "Y":
            return 24
        if chrom == "XY":
            return 25
        if chrom == "MT":
            return 26
    raise ValueError(chrom)


def default_chr_comparator(chr_a: str, chr_b: str) -> int:
    chr_a = chrom_to_int(chr_a.name)
    chr_b = chrom_to_int(chr_b.name)

    if chr_a == chr_b:
        return 0
    
    return 1 if chr_a < chr_b else -1


def fill_queue(reader: Iterable[Genotypes], queue: Queue):
    for genotypes in reader:
        queue.put(genotypes)

    queue.put(None)


def compare_variants(a: Variant, b: Variant,
                     chrom_comparator: ChromComparator = default_chr_comparator
) -> int:
    if a == b:
        return 0

    # Check if chromosome allow ordering.
    chrom_order = chrom_comparator(a.chrom, b.chrom)
    if chrom_order != 0:
        return chrom_order

    # Same chromosome, compare on positions.
    if a.pos == b.pos:
        return 0
    elif a.pos < b.pos:
        return 1
    else:
        return -1


def match_genotype_iterators(a: Iterable[Genotypes], b: Iterable[Genotypes]
) -> Iterator[Tuple[Genotypes, Genotypes]]:

    q_a = Queue()
    q_b = Queue()

    # Start processes to fill queues.
    proc_a = Process(target=fill_queue, args=(a, q_a))
    proc_b = Process(target=fill_queue, args=(b, q_b))

    processes = [proc_a, proc_b]
    for p in processes:
        p.start()

    # Main algorithm:
    # Iterate over the smallest variant.
    cur_a = q_a.get()
    cur_b = q_b.get()
    stop = False
    while not stop:
        # Iterate the last iterator until it's not last anymore.
        while compare_variants(cur_a.variant, cur_b.variant) == -1:
            # Means that B is before A.
            cur_b = q_b.get()
            if cur_b is None:
                stop = True
                break

        if stop:
            break

        while compare_variants(cur_a.variant, cur_b.variant) == 1:
            # Means that A is before B.
            cur_a = q_a.get()
            if cur_a is None:
                stop = True
                break

        if stop:
            break

        # Either variants match or same locus, different alleles.
        if cur_a.variant == cur_b.variant:
            yield cur_a, cur_b
            cur_a = q_a.get()
            cur_b = q_b.get()
        else:
            cur_a = q_a.get()  # Advance one iterator to match later.

        if cur_a is None or cur_b is None:
            stop = True

    for p in processes:
        p.join()