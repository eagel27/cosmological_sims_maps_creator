import numpy as np

phi = (1 + np.sqrt(5))/2

TETRAHEDRON_VERTICES = [(np.sqrt(8 / 9), 0, -1 / 3),
                        (-np.sqrt(2 / 9), np.sqrt(2 / 3), -1 / 3),
                        (-np.sqrt(2 / 9), -np.sqrt(2 / 3), -1 / 3),
                        (0, 0, 1)]

OCTAHEDRON_VERTICES = [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                       (0, -1, 0), (0, 0, 1), (0, 0, -1)]

CUBE_VERTICES = [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1),
                 (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1)]

ICOSAHEDRON_VERTICES = [(0, 1, phi), (0, -1, phi), (0, 1, -phi), (0, -1, -phi),
                        (1, phi, 0), (-1, phi, 0), (1, -phi, 0), (-1, -phi, 0),
                        (phi, 0, 1), (-phi, 0, 1), (phi, 0, -1), (-phi, 0, -1)]

DODECAHEDRON_VERTICES = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1),
                         (-1, -1, 1), (-1, 1, -1), (1, -1, -1), (-1, -1, -1),
                         (0, phi, 1/phi), (0, -phi, 1/phi), (0, phi, -1/phi), (0, -phi, -1/phi),
                         (1/phi, 0, phi), (-1/phi, 0, phi), (1/phi, 0, -phi), (-1/phi, 0, -phi),
                         (phi, 1/phi, 0), (-phi, 1/phi, 0), (phi, -1/phi, 0), (-phi, -1/phi, 0)]


ALIGNMENTS = {
    1: [],
    2: TETRAHEDRON_VERTICES[:2],
    3: TETRAHEDRON_VERTICES[:3],
    4: TETRAHEDRON_VERTICES,
    6: OCTAHEDRON_VERTICES,
    8: CUBE_VERTICES,
    10: ICOSAHEDRON_VERTICES[:10],
    12: ICOSAHEDRON_VERTICES,
    20: DODECAHEDRON_VERTICES
}

