#! /usr/bin/env python

from __future__ import print_function, division
from chem import collect, constants, coords 
import sys, os
from numpy import array, append, reshape
from math import sqrt

def main():
    """\
    Program for renormalizing the output of the unit sphere calculation
    in TDSPEC.

    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-i', '--initial', help="File to renormalize.", 
                        required=True)
    parser.add_argument('-r', '--renormalize', help="File to renormalize to.",
                        required=True)
    parser.add_argument('-o', '--output', help="Name of the output file.", 
                        required=True)
    parser.add_argument('-s', '--scale', help="Scale factor for the vectors.",
                        required=False,default=1.0)
    args = parser.parse_args()

    # Scale factor for the plot
    scale = float(args.scale)

    # Collect the information from the file that needs to be renormalized. 
    f = open(args.initial)
    grid = array([], dtype=float)
    contracted = array([], dtype=float)
    norm = array([], dtype=float)
    for ln in f:
        # The value used for normalizing of the file. 
        if "Maximum" in ln:
            ln = ln.rstrip().split()
            maxnorm = float(ln[4])
        # Values of the angles at the maximum.
        if "Theta(polar)" in ln:
            ln = ln.rstrip().split()
            thetaorig = float(ln[3][:-1])
            phiorig = float(ln[6])
        # The grid for plotting, and the vectors at each point. 
        if "vmd_draw_vector" in ln:
            ln = ln.rstrip().split('{')
            # The grid
            ln1 = ln[1].split('}')
            ln1 = ln1[0].split()
            # The vector 
            ln2 = ln[2].split('}')
            ln2 = ln2[0].split()
            grid = append(grid, 
                          [float(ln1[0]), float(ln1[1]), float(ln1[2])])
            contracted = append(contracted, 
                                [float(ln2[0]), float(ln2[1]), float(ln2[2])])
            norm = append(norm,
                          sqrt(maxnorm**2*(float(ln2[0])**2 +
                                           float(ln2[1])**2 +
                                           float(ln2[2])**2)))

    # Reshape the grid and contracted vectors.  There are 648 points 
    # in both arrays, because we are stepping by 10 degrees along the
    # polar angle (0 < theta < 360) and azimuthal angle (0 < phi < 180)
    # such that the number of grid points is 36 x 18 = 648.
    grid = reshape(grid, (648,3))
    contracted = reshape(contracted, (648,3))
    # Unnormalize the contracted vector.
    contracted = contracted * maxnorm

    # Collect the information from the file that we want to normalize to.
    f = open(args.renormalize)
    for ln in f:
        # The value used for normalizing of the file.
        if "Maximum" in ln:
            ln = ln.rstrip().split()
            renorm = float(ln[4])
        # Values of the angles at the maximum.
        if "Theta(polar)" in ln:
            ln = ln.rstrip().split()
            thetarenorm = float(ln[3][:-1])
            phirenorm = float(ln[6])

    # Renormalize the contracted vector.
    contracted = contracted / renorm 

    # Find the new maximum due to renormalization (this won't be 1, in
    # general).
    newmax = maxnorm / renorm

    # Output some data to the new file.
    f = open(args.output, 'w')
    f.write('# RENORMALIZED FILE\n')
    f.write('# Original normalization information\n')
    f.write('{0} {1:12.5E}{2}'.format('# Original Normalization Maximum (a.u.):', maxnorm, '\n'))
    f.write('{0} {1:5.1F}{2} {3:5.1F}{4}'.format('# Theta(polar) =', thetaorig, 
                                                 ', Phi(azimuthal) =', phiorig, '\n'))
    f.write('# New normalization information\n')
    f.write('{0} {1:12.5E}{2}'.format('# New Normalization Maximum (a.u.):', renorm, '\n'))
    f.write('{0} {1:5.1F}{2} {3:5.1F}{4}'.format('# Theta(polar) =', thetarenorm, 
                                                 ', Phi(azimuthal) =', phirenorm, '\n'))
    f.write('{0} {1:8.5F}{2}'.format('# Original Maximum to New Maximum Ratio =', newmax, '\n'))
    f.write('# This shows the new maximum resulting from renormalization\n')
    f.write('{0} {1:6.4F}{2}'.format('# Scale factor for plot:', scale, '\n'))
    f.write('color change rgb 11 0.1 0.1 1.0\n')
    f.write('color change rgb 12 0.2 0.2 1.0\n') 
    f.write('color change rgb 13 0.3 0.3 1.0\n')
    f.write('color change rgb 14 0.4 0.4 1.0\n')
    f.write('color change rgb 15 0.5 0.5 1.0\n')
    f.write('color change rgb 17 0.6 0.6 1.0\n')
    f.write('color change rgb 18 0.7 0.7 1.0\n')
    f.write('color change rgb 19 0.8 0.8 1.0\n')
    f.write('color change rgb 20 0.9 0.9 1.0\n')
    f.write('color change rgb 21 1.0 0.1 0.1\n')
    f.write('color change rgb 22 1.0 0.2 0.2\n')
    f.write('color change rgb 23 1.0 0.3 0.3\n')
    f.write('color change rgb 24 1.0 0.4 0.4\n')
    f.write('color change rgb 25 1.0 0.5 0.5\n')
    f.write('color change rgb 26 1.0 0.6 0.6\n')
    f.write('color change rgb 27 1.0 0.7 0.7\n')
    f.write('color change rgb 28 1.0 0.8 0.8\n')
    f.write('color change rgb 29 1.0 0.9 0.9\n')
    fmt = '{0}{1:8.5F} {2:8.5F} {3:8.5F}{4}{5:12.5E} {6:12.5E} {7:12.5E}{8} {9:5.3F}{10}'
    for i in range(648):
        color = norm[i] / renorm
        if (0.00000<=color and color<0.047619*newmax):
            f.write('draw color 0\n') # Blue (smallest intensity)
        elif (0.047619*newmax<=color and color<0.095238*newmax):
            f.write('draw color 11\n')
        elif (0.095238*newmax<=color and color<0.142857*newmax):
            f.write('draw color 12\n')
        elif (0.142857*newmax<=color and color<0.190476*newmax):
            f.write('draw color 13\n')
        elif (0.190476*newmax<=color and color<0.238095*newmax):
            f.write('draw color 14\n')
        elif (0.238095*newmax<=color and color<0.285714*newmax):
            f.write('draw color 15\n')
        elif (0.285714*newmax<=color and color<0.333333*newmax):
            f.write('draw color 17\n')
        elif (0.333333*newmax<=color and color<0.380952*newmax):
            f.write('draw color 18\n')
        elif (0.380952*newmax<=color and color<0.428571*newmax):
            f.write('draw color 19\n')
        elif (0.428571*newmax<=color and color<0.476190*newmax):
            f.write('draw color 20\n')
        elif (0.476190*newmax<=color and color<0.523809*newmax):
            f.write('draw color 8\n') # White
        elif (0.523809*newmax<=color and color<0.571428*newmax):
            f.write('draw color 29\n')
        elif (0.571428*newmax<=color and color<0.619047*newmax):
            f.write('draw color 28\n')
        elif (0.619047*newmax<=color and color<0.666666*newmax):
            f.write('draw color 27\n')
        elif (0.666666*newmax<=color and color<0.714285*newmax):
            f.write('draw color 26\n')
        elif (0.714285*newmax<=color and color<0.761904*newmax):
            f.write('draw color 25\n')
        elif (0.761904*newmax<=color and color<0.809523*newmax):
            f.write('draw color 24\n')
        elif (0.809523*newmax<=color and color<0.857142*newmax):
            f.write('draw color 23\n')
        elif (0.857142*newmax<=color and color<0.904761*newmax):
            f.write('draw color 22\n')
        elif (0.904761*newmax<=color and color<0.952380*newmax):
            f.write('draw color 21\n')
        elif (0.952380*newmax<=color and color<=1.000000*newmax):
            f.write('draw color 1\n') # Red (largest intensity)
        f.write(fmt.format('vmd_draw_vector 0 {', grid[i][0],
                           grid[i][1], grid[i][2], '} {', 
                           contracted[i][0], contracted[i][1],
                           contracted[i][2], '}', scale, ' 30 0.08\n'))
    f.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
