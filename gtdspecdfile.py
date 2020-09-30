#! /usr/bin/env python

from __future__ import print_function, division
import sys, os

def main():
    """\
    This script projects the normal modes onto either the excited state
    optimized geometry or excited state gradient.  The quantities required
    are evaluated analytically in ADF v2010 or later, allowing for hasty
    generation of deltas for resonance Raman and absorption spectra
    simulations.

    The resulting dimensionless deltas are written to screen or file in the
    for of a TDSPEC datafile, with all fields filled based on data from the
    collected file.  Note that since the 'gamma' parameter is phenomenolgical
    a default of 400 will be used, which may be changed at the users
    discretion with the gamma command line options.

    Only two files are required to calculate the deltas:

    1. ADF output file for a frequencies calculation, done either numerically
       or analytically.
    2. ADF output file containing the analytical excited state gradients or
       optimized excited state coordinates.

    Two methods exist that allow you to evaluate the dimensionless 
    displacements using analytical gradient methods:

    1. Projecting normal modes onto the difference between the ground and
    excited state equilibrium geometry.  This will be called the analytical 
    excited state optimization method.

    2. Projecting the normal modes onto the analytical gradient of the 
    excited state energy at the ground state equilibrium position.  This
    will be called the analytical excited state gradient method.

    This program has the capability to use both methods, and can be used 
    to test both methods against each other.  In terms of calculation 
    speed, the latter is faster because you do not need the analytical 
    geometry optimization calculation to finish, since the first step of 
    the geometry optimization gives the required gradient.

    Note that the program determines what delta to calculate based on if
    a full excited state optimization or just the gradient was calculated
    (based on if the number of geometry iterations is greater than one or
    not).  You may explicity determine which one to use with a command line
    option.

    Additionally, TDSPEC may use more than one excited state.  Therefore,
    you may give more than one ES file to make a datafile containing multiple
    excited states.

    AUTHORS:

    Daniel Silverstein
    Seth M. Morton

    CHANGELOG:

    6-29-2012 DWS v1.1

    Added scaling of vibrational frequencies and derivatives of the ground
    state dipole moment and transition dipole moment for generation of 
    various kinds of spectra.

    2-1-2011 SMM v1.0

    Migrated from perl to python.  Reduced required input files from
    4 to 2.
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from argparse import FileType
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('modes', help='The output file for a frequencies '
                        'calculation.')
    parser.add_argument('excited', help='The excited state output file(s) to '
                        'use to calculate the dimensionless deltas.  They may '
                        'be either full optimizations or just the gradients. '
                        'Note that the ground state coordinates are also '
                        'taken from this file.  If doing numerical deltas, '
                        'these are just the equilibrium excitations, and '
                        'therefore only one of these files is required (even '
                        'for multiple deltas) since they all share the same '
                        'equilibrium position.', nargs='+')
    parser.add_argument('--output', '-o', help='Name of the file to write to. '
                        ' If none is given, the datafile will be written to '
                        'screen.', type=FileType('w'), default=sys.stdout)
    parser.add_argument('-g', '--gamma', help='Set the value of the '
                        'phenomenological damping parameter "gamma".  The '
                        'default is %(default)f.', default=400.0, type=float)
    parser.add_argument('--esgrad', help='Force this program to calculate '
                        'deltas using the excited state gradients, not the '
                        'excited state optimized geometry.', default=False,
                        action='store_true')
    parser.add_argument('--excite', '-e', help='The excitations to examine '
                        "for numerical deltas.  Must be in '# sym' form, "
                        'including quotes.', nargs='+')
    parser.add_argument('--dirs', '-d', help='The directories where the '
                        'numerical delta data will be taken from.  If you '
                        'are calculating several numerical deltas and they '
                        'all come from the same source, you only need to '
                        'list the one source, and it will be used for all '
                        'deltas', nargs='+')
    parser.add_argument('--scale', '-s', help='Scale the vibrational '
                        'frequencies by a factor.', type=float, default=1.0)
    parser.add_argument('--IR', '-i', help='Calculate derivatives of the '
                        'ground state dipole moment for an SFG spectrum.', 
                        default=False, action='store_true')
    parser.add_argument('--dtdip', '-t', help='Calculate derivatives of the '
                        'transition dipole moment for HT terms.', default=False,
                        action='store_true')
    parser.add_argument('--dstpm', '-S', help='Calculate derivatives of the '
                        'two-photon transition moment for RHRS B2 terms.', 
                        default=False, action='store_true')
    parser.add_argument('--low', help='The low mode to include. The '
                        'default is %(default)s cm-1.', default=float(400),
                        type=float)
    parser.add_argument('--high', help='The high mode to include. The '
                        'default is %(default)s cm-1.', default=float(2000),
                        type=float)
    args = parser.parse_args()

    #print(args.IR, args.dtdip, args.dstpm)

    # Collect information from the given files (modes contains the
    # normal mode file, es is the excitations calculation for the 
    # optimized geometry, and also has the dimensionless displacements.
    modes, es = collect_values(args)

    # Calculate the deltas, depending on the method requested by the
    # user. Xing
    modes = calculate_deltas(modes, es, args)

    # Print the datafile xing
    print_datafile(modes, es, args)

def collect_values(args):
    '''Collect the pertinent information from file, and verify it.'''
    from chemPackage import collect
    from numpy import array, reshape

    # Check that when the excited states are listed with the excite option,
    # that we also list a directory for collecting data from. 
    if (args.excite and not args.dirs) or (args.dirs and not args.excite):
        sys.exit('Need either both or neither of --excite and --dirs')
    # Verify that we have the correct number of directories to collect
    # data from.
    if args.excite:
        if len(args.excite) != len(args.dirs) and  len(args.dirs) > 1:
            sys.exit('# of --excite and --dirs must be equal')
        # If --dirs is 1, then fill to match excite
        if len(args.dirs) == 1 and len(args.excite) > 1:
            for i in range(1, len(args.excite), 1):
                args.dirs.append(args.dirs[0])
        # If --excited is 1, then fill to match excite
        if len(args.excited) == 1 and len(args.excite) > 1:
            for i in range(1, len(args.excite), 1):
                args.excited.append(args.excited[0])
            
    # Collect the ground state geometry and normal modes, and make
    # sure it is of the correct type
    modes = collect(args.modes)
    assert 'FREQUENCIES' in modes.calctype, (
       'The MODES file {0} has no frequencies info.'.format(modes.filename))

    # Collect all the excited state information, and make
    # sure it is of the correct type
    es = []
    for i, f in enumerate(args.excited):
        if args.excite:
            # Collect the equilibrium excitation
            es.append(collect(f))
            # Determine the sign of the transition dipole moment and two-photon
            # transition moment at the equilibrium geometry
            sign_tdm, sign_stpm = transmom_sign(es[i], args)
            # Collect the numerical deltas from the given directories and the
            # given excitations.
            temp = modes.copy()
            temp.collect_raman_derivatives(dir=args.dirs[i],
                                           excitation=args.excite[i],
                                           gdipder=args.IR,
                                           tdipder=args.dtdip,
                                           tdipref=sign_tdm)
            #if args.dstpm:
            #
            # Place the delta and excited state number information into that dataset
            es[i].include(temp)
            es[i].es = args.excite[i]
            assert 'DELTAS' in es[i].calctype, (
                    'The EXCITED directory {0} did not yield deltas'.format(
                                                   args.dirs[i]))
        else:
            es.append(collect(f))
            assert 'EXCITED STATE' in es[i].calctype, (
                               'The EXCITED file {0} has no es info.'.format(f))

    # Make the number of normal modes match the number of numerical deltas
    if args.excite:
        # If the number of deltas is not the same for all, exit
        if len(set([x.nmodes for x in es])) != 1:
            sys.exit('Each set of numerical deltas must be of the same length')
        else:
            modes.nmodes        = es[0].nmodes
            modes.IR            = es[0].IR.copy()
            modes.normal_modes  = es[0].normal_modes.copy()
            modes.v_frequencies = es[0].v_frequencies.copy()

    # Determine if we will use the excited state gradient or optimized
    # geometry method.  Base this on the first file given, and only if
    # not chosen explicitly on the command line.
    if not args.esgrad:
        args.esgrad = True if 'GEOMETRY' not in es[0].calctype else False
   
    return modes, es

def calculate_deltas(md, es, args):
    '''Prep for the Delta calculation, and determine which method to use.'''
    from chemPackage.constants import atomic_mass, HBAR, PI, LIGHT
    from numpy import fastCopyAndTranspose as fcat
    from numpy import array, sqrt, zeros

    # Define the conversion factor
    CONVERSION = sqrt(HBAR / ( 2 * PI * LIGHT * 100 ))

    # Create an array that holds all the root masses of each atom
    rt_mass = fcat(array([sqrt(md.masses), sqrt(md.masses), sqrt(md.masses)]))

    # Allocate the deltas
    md.deltas = zeros((len(es), md.nmodes), dtype=float)
    md.dgdip = zeros((len(es),md.nmodes,3), dtype=float)
    if args.dtdip:
        for n, e in enumerate(es):
            # Because we sometimes use a different number of excited states
            # for the equilibrium geometry calculation than the gradient
            # calculations, we need to check that the Numpy arrays have the
            # same dimension.  This prevents a lot of errors.
            dim = e.dtdip.shape[1]
            if dim != e.nexcite:
                md.dtdip = zeros((len(es),md.nmodes,dim,3), dtype=float)
            else:
                md.dtdip = zeros((len(es),md.nmodes,e.nexcite,3), dtype=float)

    # Choose the method to calculate/collect the deltas
    # We also scale the frequencies.
    md.scale_vfreq(args.scale)
    if args.excite:
        # Numerical gradients.
        for n, e in enumerate(es): 
            md.deltas[n] = e.deltas.copy()
            if args.IR:
                md.dgdip[n] = e.dgdip.copy()        
            if args.dtdip:
                md.dtdip[n] = e.dtdip.copy()
#        md.scale_vfreq(args.scale)
        return md
    elif args.esgrad:
        # Analytical excited state gradients.
        return gradient_method(md, es, CONVERSION, rt_mass, args)
    else:
        # Deltas based on excited state geometry optimization.
        return optimization_method(md, es, CONVERSION, rt_mass, args)

def optimization_method(md, es, CONVERSION, rt_mass, args):
    '''\
    Based on the review:

    F. Neese, T. Petrenko, D. Ganyushin, and G. Olbrich.  Coord. Chem. 
    Rev., 251, 288 (2010).

    The analytical excited state optimization method is applied by using
    the equation (see page 318, Eqs. 108 and 109 of the review):

    dim_delta_j = sum_i[L_{ij}*Dm_i]*SQRT(omega_j/hbar)

    Where L_{ij} is the normal mode vector in normal coordinates, Dm_i is 
    the mass-weighted difference between the GS and ES geometry, and 
    omega_j is the vibrational frequency (omega_j = 2*PI*c*nu_j).

    **NOTE** the mass-weighted delta has units of mass^(1/2)*length.
    '''
    from chemPackage.constants import AMU
    from numpy import sum as npsum
    from numpy import zeros
    from math import sqrt

    # Loop over each excited state calculation
    for n in range(len(es)):

        # Find the mass-weighted difference between the GS and ES geometry,
        # then convert to SI units
        Dm = (es[n].coordinates - es[n].initial_coordinates) * rt_mass
        Dm *= sqrt(AMU) * 1E-10

        # For each normal mode, perform the equation 
        # dim_delta_j = sum_i[L_{ij}*Dm_i]*SQRT(omega_j/hbar)
        # Note that i is each atom, but this loop is done implicitly with npsum
        for j in range(md.nmodes):

            # Convert the normal modes into normal coordinates (L_{ij}).
            # We use the absolute value in case of imaginary frequencies
            try:
                sqrt_omega_hbar = CONVERSION / sqrt(abs(md.v_frequencies[j]))
                md.normal_modes[j] *= sqrt(AMU) * 1E-10 * rt_mass / sqrt_omega_hbar

                # Now actually calculate the dimensionless deltas
                md.deltas[n][j] += npsum(md.normal_modes[j] * Dm / sqrt_omega_hbar)
            except FloatingPointError:
                md.deltas[n][j] += 0.

        # Scale the vibrational frequency
#        md.scale_vfreq(args.scale)

    return md

def gradient_method(md, es, CONVERSION, rt_mass, args):
    '''\
    Based on the review:

    F. Neese, T. Petrenko, D. Ganyushin, and G. Olbrich.  Coord. Chem. Rev., 
    251, 288 (2010).

    The analytical excited state gradient method is applied by using the equation
    (see p. 319, Eqs. 113 and 114 in the review):

    dim_delta_j = -1/omega_j*SQRT(hbar/omega_j)*V_{Q_j}
    V_{Q_j} = sum_i[(1/sqrt(m_i))*V_{X_i}*L_{ij}]

    Where omega_j is the vibrational frequency (omega_j = 2*pi*c*nu_j), V_{Q_j}
    is the excitation energy gradient along mode Q_j in mass-weighted units,
    and L_{ij} is the normal mode vector in normal coordinates.
    '''
    from chemPackage.constants import AMU, PI
    from chemPackage.constants import HART2WAVENUM as H2WN
    from chemPackage.constants import ANGSTROM2BOHR as A2B
    from chemPackage.constants import BOHR2CM as B2CM
    from numpy import sum as npsum
    from numpy import zeros
    from math import sqrt

    # Loop over each excited state calculation
    for n in range(len(es)):

        # Convert Cartesian energy gradient into mass-weighted energy gradient
        # Units here are [Eh/(Angstrom * sqrt(Dalton))] 
        grad = es[n].es_gradient['total TDDFT gradient'] * A2B / rt_mass

        # Pseudo-atomic units (mass is reported in amu) [Eh/(a0*sqrt(Dalton))]
        #grad = es[n].es_gradient['total TDDFT gradient'] / rt_mass

        for j in range(md.nmodes):

            # Convert the normal modes into normal coordinates (L_{ij}).
            # We use the absolute value in case of imaginary frequencies
            # The units for the normal mode are dimensionless.
            try:
                sqrt_omega_hbar = CONVERSION / sqrt(abs(md.v_frequencies[j]))
                md.normal_modes[j] *= sqrt(AMU) * 1E-10 * rt_mass / sqrt_omega_hbar

                mwgrad = npsum(grad * md.normal_modes[j])

                # Calculate the dimensionless deltas
                # Extra factor of 2*pi in the fourth term comes from the fact that
                # this is an angular frequency.  We need the vibrational frequency
                # to cancel out energy units, so the speed of light is not
                # multiplied in.  The absolute value is used to avoid errors due to
                # imaginary frequencies.
                md.deltas[n][j] = ( - mwgrad
                                  * ( 1 / sqrt(AMU) )
                                  * ( 1 / 1E-10 )
                                  * ( H2WN / ( 2 * PI * md.v_frequencies[j] ) )
                                  * ( CONVERSION / sqrt(abs(md.v_frequencies[j])) )
                                  )
            except FloatingPointError:
                md.deltas[n][j] = 0.

            ## Neese and Petrenko method
            ## Units of the normal mode vector are [a0*sqrt(Dalton)]
            #md.normal_modes[j] *= rt_mass * A2B

            ## Units of the mass-weighted displacements are [hbar^2 / Eh]
            ## These units are mass*length^2 
            #mwdelta = npsum(grad * md.normal_modes[j])

            #if (md.v_frequencies[j] > 0.0):
            #    mwdelta *= ( H2WN * H2WN )/( md.v_frequencies[j] * md.v_frequencies[j] )
            #else:
            #    mwdelta = 0.0

            #mwdelta = sqrt(abs(mwdelta))

            ##md.deltas[n][j] = ( mwdelta
            ##                  * ( 1 / 0.000548579903 )
            ##                  * (( md.v_frequencies[j] * md.v_frequencies[j] )/( conv * conv ))**0.25 )
            ##                  * md.v_frequencies[j] * md.v_frequencies[j] * B2CM * B2CM )

            ## Convert to SI units [kg*m^2]
            ##mwdelta *= B2CM * B2CM * 0.01 * 0.01 * AMU             
            #mwdelta *= B2CM * 0.01 * sqrt(AMU)

            #md.deltas[n][j] = ( mwdelta 
            #                  * ( sqrt(abs(md.v_frequencies[j])) / CONVERSION )
            #                  )

        # Scale the vibrational frequency
#        md.scale_vfreq(args.scale)

    return md

def transmom_sign(es, args):
    '''\
    Determine the sign of the largest component of the transition
    moment.
    '''
    from numpy import absolute, argmax

    sign_tdm = []
    sign_stpm = []
    
    # For situations when both the derivative of the transition dipole
    # moment and two-photon transition moment are needed.
    if args.dtdip and args.dstpm:
        temp1 = absolute(es.TDM)
        temp1 = temp1.argmax(axis=1)
        temp2 = absolute(es.STPM)
        temp2 = temp2.argmax(axis=1)
    elif args.dtdip:
        temp1 = absolute(es.TDM)
        temp1 = temp1.argmax(axis=1)

    for num in range(es.nexcite):    
        if args.dtdip and args.dstpm:
            sign_tdm.append(es.TDM[num][temp1[num]]/absolute(es.TDM[num][temp1[num]]))
            sign_stpm.append(es.STPM[num][temp2[num]]/absolute(es.STPM[num][temp2[num]]))
        elif args.dtdip:
            sign_tdm.append(es.TDM[num][temp1[num]]/absolute(es.TDM[num][temp1[num]]))

    return sign_tdm, sign_stpm
    

def print_datafile(md, es, args):
    '''Print the datafile to file or screen.'''
    from numpy import where

    # Format strings
    d  = '{0:>7.2f} {1:11.8f}'
    sfgarrsb = '{0:>7.2f} {1:11.8f} {2:8.5f} {3:8.5f} {4:8.5f}'
    sfgb = '{0:>7.2f} {1:11.8f} {2:8.5f} {3:8.5f} {4:8.5f} {5:8.5f} {6:8.5f} {7:8.5f}'
    t  = 'Tdip {0[0]:8.5f} {0[1]:8.5f} {0[2]:8.5f}'
    en = 'Energy {0:10.8f}'

    # The total number of excitations
    print('NEXCI {0:d}'.format(len(es)), file=args.output)

    # For each excitation, print parameters and deltas
    for n, e in enumerate(es):

        # This bit is a temp fix for the dtdip for higher excitations
        if n == 0:
            iex = args.excite

        # Energy of the excitation.  Find the excitation that matches what was
        # selected in the input.
#        num, sym = e.es.split()
#        indx = where(e.excitation_symmetries == sym)[0]
#        num = int(num) - 1
#        if args.excite:
        try:
            indx = where(e.excitation_symmetries == e.excite)[0]
            es_energy = e.excitation_energies[indx]
            tdip = e.TDM[indx]
        except AttributeError:
#        else:
            num, sym = e.es.split()
            num = int(num) - 1
            indx = where(e.excitation_symmetries == sym)[0][num]
            es_energy = e.excitation_energies[indx]
            tdip = e.TDM[indx]

        # The number of this excitation
        if args.excite:
            print('Excitation {0:d}'.format(n+1), file=args.output)
        else:
            print('Excitation {0:d}'.format(num+1), file=args.output)

        print(en.format(es_energy), file=args.output)

        # Damping parameter
        print('Gamma {0:6.1f}'.format(args.gamma), file=args.output)

        # Transition dipole moment of the excitation
        print(t.format(tdip), file=args.output)

        # The number of frequencies
        nfreq = 0
        for j in range(md.nmodes):
            freq = md.v_frequencies[j]
            if (freq > args.low) and (freq < args.high):            
                nfreq += 1
        #print('Frequencies {0:d}'.format(md.nmodes), file=args.output)
        print('Frequencies {0:d}'.format(nfreq), file=args.output)
   
        # Print off each delta for each mode of this excitation.
        for j in range(md.nmodes):
            freq = md.v_frequencies[j]
            if (freq > args.low) and (freq < args.high):            
                if (args.IR and args.dtdip):
                    print(sfgb.format(md.v_frequencies[j], 
                                          md.deltas[n][j],
                                          e.dgdip[j][0],
                                          e.dgdip[j][1],
                                          e.dgdip[j][2],
                                          e.dtdip[j][n][0],
                                          e.dtdip[j][n][1],
                                          e.dtdip[j][n][2]),
                          file=args.output)
                elif args.IR:
                    print(sfgarrsb.format(md.v_frequencies[j],
                                          md.deltas[n][j],
                                          e.dgdip[j][0],
                                          e.dgdip[j][1],
                                          e.dgdip[j][2]),
                          file=args.output)
                elif args.dtdip:
                    print(sfgarrsb.format(md.v_frequencies[j], 
                                          md.deltas[n][j],
                                          e.dtdip[j][n][0],
                                          e.dtdip[j][n][1],
                                          e.dtdip[j][n][2]),
                          file=args.output)
                else:
                    print(d.format(md.v_frequencies[j], md.deltas[n][j]),
                          file=args.output)

        # Print the end of the block
        print('End', file=args.output)


if __name__ == '__main__':
    try:
        main()
    except AssertionError as a:
        sys.exit(a)
    except KeyboardInterrupt:
        sys.exit(1)
