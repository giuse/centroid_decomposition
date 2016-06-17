# Memory-Efficient Centroid Decomposition

As seen in:

- Mourad Khayati, Michael H. Böhlen, and Johann Gamper. “Memory-Efficient Centroid Decomposition for Long Time Series.” In IEEE 30th International Conference on Data Engineering, Chicago, ICDE 2014, IL, USA, March 31 - April 4, 2014, 100–111, 2014. Bibtex PDF

More info:

- http://exascale.info/members/mourad-khayati/

## Paper errata:

- Pag 105, algo SSV, search nexe element:
    `vi.abs > val` should be `vi.abs > val.abs`
- Pag 106, example 3, second iteration:
    use Z(2) to compute **_both S(2) and_** V(2)
