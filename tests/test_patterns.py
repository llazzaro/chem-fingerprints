
hydrogen_test_cases = [
    ("[C]", 0),
    ("[CH]", 1),
    ("[CH2]", 2),
    ("[CH3]", 3),
    ("[CH4]", 4),
    ("[C][C]", 0),
    ("[C][CH]", 1),
    ("[CH][CH2]", 3),
    ("[H]", 1),
    ("[HH]", 2),
    ("[H][CH3]", 4),
    ("C", 4),
    ("[CH]([2H])([3H])[H]", 4),
    ("[H][CH3]", 4),
    ("[H][H]", 2),
    ]

aromatic_ring_cases = [
    ("C1CCCCC1", 0),
    ("c1ccccc1", 1),
    ("c1ccccc1.c1ccccc1", 2),
    ("c1cncc2c1nncn2", 2),
    ("c1csc2c1csc2", 2),
    ("c1csc2c1csc2.c1csc2c1csc2", 4),
    ("c1ccccc1.c1ccccc1.c1ccccc1", 3),
    ("c1ccccc1.c1ccccc1.c1ccccc1.c1ccccc1", 4),
    ("c1ccccc1.c1ccccc1.c1ccccc1.c1ccccc1.c1ccccc1", 5),
    ("c1ccccc1.c1ccccc1.c1ccccc1.c1ccccc1.c1ccccc1.c1ccccc1", 6),
    ]
