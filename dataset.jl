using PolaronMobility

MAPI = material(4.9, 24.1, 0.12, 2.25)

p = polaron(MAPI, verbose = true)

p_2 = polaron(MAPI, [1, 100, 200], [1,2,3,4], verbose = true)