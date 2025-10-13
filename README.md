# Motif_Maximization_Analysis-
This directory represents a maximization analysis where the goal is to take a set of SCRMshaw Cis-regulatory modules (CRMs),
and find the common motifs seen in these regions as well as a maximization/co-occurrence analysis

The analysis works by running the CRM set through FIMO, receiving a fimo.tsv file which includes all possible motifs.
Through fimo.tsv we can begin studying co-occurrence. By implementing a motif counting script (Motif_Counter_per_seq.py),
on the fimo.tsv file, we receive a dictionary where the key is the Motif ID, and the values are the number of times the motif
appears in the total set, and the number of sequences the motifs appear in that are in the set.

Using this dictionary we can first find the most commonly found motif in the CRM set. We then proceed to use all of the sequences,
that the most common motif, Motif 1, is found in, as the input for a second FIMO run. This second FIMO run will generate a fimo.tsv,
where we are sure that Motif 1 is in, but now we can find the 2nd most commonly occurring motif in this set, Motif 2. This pattern
is repeated an arbitrary number of times. The goal is to come down to a set of motifs that are a commonly found together. Motif 2 is
found in a lot of or all sequences of motif 1. Motif 3 is found in a lot of or all sequences of Motif 2 AND Motif 1. This creates
a set of sequences that we now know have this list of commonly occurring motifs, [Motif 1, Motif 2, Motif 3]

Questions to ask:

Do different sets of CRMs from different samples, have similarly occurring motifs?

Do different sets of CRMs from different samples, have a similar ordering of motifs?

How common are these motif occurrences within each CRM set?

How unique are motif occurrences between each CRM set?


A minimum number of scripts/files needed for this analysis:

        1. A motif database
        2. A bash fimo script to run FIMO
        3. A motif counting script
