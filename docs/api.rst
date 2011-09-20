.. _chemfp-api:

==========
chemfp API
==========

This chapter contains the docstrings for the public portion of the
chemfp API.

(Note: the values of 0.69999999999999996 which you see are actually
0.7. I don't know how to get my documentation tool to use the more
understandable value.)

.. py:module:: chemfp

.. autofunction:: chemfp.open
.. autofunction:: chemfp.load_fingerprints

.. _chemfp_read_structure_fingerprints:
.. autofunction:: chemfp.read_structure_fingerprints

.. _chemfp_count_tanimoto_hits:
.. autofunction:: chemfp.count_tanimoto_hits

.. _chemfp_threshold_tanimoto_search:
.. autofunction:: chemfp.threshold_tanimoto_search
.. autofunction:: chemfp.knearest_tanimoto_search

.. autoclass:: chemfp.Metadata
.. autoclass:: chemfp.FingerprintIterator
.. autoclass:: chemfp.Fingerprints

.. _FingerprintArena:

.. autoclass:: chemfp.arena.FingerprintArena
.. autoclass:: chemfp.arena.SearchResults

.. _chemfp.bitops:
.. py:module:: chemfp.bitops

.. autofunction:: chemfp.bitops.byte_popcount
.. autofunction:: chemfp.bitops.byte_intersect_popcount
.. autofunction:: chemfp.bitops.byte_tanimoto
.. autofunction:: chemfp.bitops.byte_contains
.. autofunction:: chemfp.bitops.hex_isvalid
.. autofunction:: chemfp.bitops.hex_popcount
.. autofunction:: chemfp.bitops.hex_intersect_popcount
.. autofunction:: chemfp.bitops.hex_tanimoto
.. autofunction:: chemfp.bitops.hex_contains
