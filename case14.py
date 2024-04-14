# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 18:04:09 2023

@author: ULRICH_LUMENDO
"""  
from numpy import array
def case():
    ppc = {"version": '2'}
    ppc["baseMVA"] = 100.0
    ppc["bus"] = array([
    		[1.0,    3.0, 0.0,      0.0,      0.0,    0.0,    1.0, 1.06,     0.0,      0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[2.0,    2.0, 21.7,     12.7,     0.0,    0.0,    1.0, 1.045,    -4.98,    0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[3.0,    2.0, 94.2,     19.0,     0.0,    0.0,    1.0, 1.01,     -12.72,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[4.0,    1.0, 47.8,     -3.9,     0.0,    0.0,    1.0, 1.019,    -10.33,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[5.0,    1.0, 7.6,      1.6,      0.0,    0.0,    1.0, 1.02,     -8.78,    0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[6.0,    2.0, 11.2,     7.5,      0.0,    0.0,    1.0, 1.07,     -14.22,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[7.0,    1.0, 0.0,      0.0,      0.0,    0.0,    1.0, 1.062,    -13.37,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[8.0,    2.0, 0.0,      0.0,      0.0,    0.0,    1.0, 1.09,     -13.36,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[9.0,    1.0, 29.5,     16.6,     0.0,    19.0,   1.0, 1.056,    -14.94,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[10.0,   1.0, 9.0,      5.8,      0.0,    0.0,    1.0, 1.051,    -15.1,    0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[11.0,   1.0, 3.5,      1.8,      0.0,    0.0,    1.0, 1.057,    -14.79,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[12.0,   1.0, 6.1,      1.6,      0.0,    0.0,    1.0, 1.055,    -15.07,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[13.0,   1.0, 13.5,     5.8,      0.0,    0.0,    1.0, 1.05,     -15.16,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ],
    		[14.0,   1.0, 14.9,     5.0,      0.0,    0.0,    1.0, 1.036,    -16.04,   0.0,      1.0, 1.06,     0.94,     0.6, 10   ]])
    ppc["gen"] = array([
    		[1.0,    232.4,    -16.9,    10.0,     0.0,      1.06,     100.0,  1.0, 332.4,    0.0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 99.72,    132.96    ],
    		[2.0,    40.0,     42.4,     50.0,     -40.0,    1.045,    100.0,  1.0, 140.0,    0.0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 42.0,     56.0      ],
    		[3.0,    0.0,      23.4,     40.0,     0.0,      1.01,     100.0,  1.0, 100.0,    0.0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30.0,     40.0      ],
    		[6.0,    0.0,      12.2,     24.0,     -6.0,     1.07,     100.0,  1.0, 100.0,    0.0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30.0,     40.0      ],
    		[8.0,    0.0,      17.4,     24.0,     -6.0,     1.09,     100.0,  1.0, 100.0,    0.0,      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30.0,     40.0      ]])
    ppc["branch"] = array([
    		[1.0,    2.0,    0.01938,       0.05917,       0.0528,        9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[1.0,    5.0,    0.05403,       0.22304,       0.0492,        9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[2.0,    3.0,    0.04699,       0.19797,       0.0438,        9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[2.0,    4.0,    0.05811,       0.17632,       0.034,         9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[2.0,    5.0,    0.05695,       0.17388,       0.0346,        9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[3.0,    4.0,    0.06701,       0.17103,       0.0128,        9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[4.0,    5.0,    0.01335,       0.04211,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[4.0,    7.0,    0.0,           0.20912,       0.0,           9900.0, 0.0,0.0,0.978,  0.0,1.0,-360.0, 360.0,  0.1     ],
    		[4.0,    9.0,    0.0,           0.55618,       0.0,           9900.0, 0.0,0.0,0.969,  0.0,1.0,-360.0, 360.0,  0.1     ],
    		[5.0,    6.0,    0.0,           0.25202,       0.0,           9900.0, 0.0,0.0,0.932,  0.0,1.0,-360.0, 360.0,  0.1     ],
    		[6.0,    11.0,   0.09498,       0.1989,        0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[6.0,    12.0,   0.12291,       0.25581,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[6.0,    13.0,   0.06615,       0.13027,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[7.0,    8.0,    0.0,           0.17615,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[7.0,    9.0,    0.0,           0.11001,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[9.0,    10.0,   0.03181,       0.0845,        0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[9.0,    14.0,   0.12711,       0.27038,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[10.0,   11.0,   0.08205,       0.19207,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[12.0,   13.0,   0.22092,       0.19988,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ],
    		[13.0,   14.0,   0.17093,       0.34802,       0.0,           9900.0, 0.0,0.0,0.0,    0.0,1.0,-360.0, 360.0,  0.1     ]])
    ppc["gencost"] = array([
    		[2.0, 0.0, 0.0, 3.0, 0.0430293,20.0,   0.0,    24.0,   12.0    ],
    		[2.0, 0.0, 0.0, 3.0, 0.25,20.0,   0.0,    24.0,   12.0    ],
    		[2.0, 0.0, 0.0, 3.0, 0.01,40.0,   0.0,    48.0,   24.0    ],
    		[2.0, 0.0, 0.0, 3.0, 0.01,40.0,   0.0,    48.0,   24.0    ],
    		[2.0, 0.0, 0.0, 3.0, 0.01,40.0,   0.0,    48.0,   24.0    ]])
    return ppc