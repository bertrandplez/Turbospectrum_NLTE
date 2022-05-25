#!/bin/csh -f

#
# This script is a demo script for NLTE TS.
#

date
set mpath = ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/COM/TEST-data
set dpath = ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/COM/TEST-data

set MODEL = atmos.sun_marcs_t5777_4.44_0.00_vmic1_new

set Caabu = 6.34
set Feabu = 7.50

set lam_min    = '4800.'
set lam_max    = '6800.'

set deltalam   = '0.006'
set METALLIC   = '     0.000'
set TURBVEL = 1.0

time ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/exec/babsma_lu <<EOF
###########
# wavelength range for the continuous opacity calculations. Should encompass the full 
# range asked for in the following spectrum calculation (bsyn_lu)
# the step is set to 1A in babsma if smaller than 1A here.
#
'LAMBDA_MIN:'  '${lam_min}'
'LAMBDA_MAX:'  '${lam_max}'
'LAMBDA_STEP:' '${deltalam}'
###########
# model atmosphere. Various formats allowed. Only MARCS can be binary or ascii. Others are ascii
#
'MODELINPUT:' '$mpath/${MODEL}'
'MARCS-FILE:' '.false.'
###########
# output continuous opacity file providing continuous abs and scatt at all atmospheric depths
# for a set of wavelengths defined by lambda_min/max/step. If the step is < 1A, it is set to 
# 1A by default 
#
'MODELOPAC:' 'contopac/${MODEL}opac'
###########
# Chemical composition. First overall metallicity, then alpha/Fe, Helium/H, and r- and s-process
# The latter are scaled according to their solar-system fraction (see makeabund.f)
# finally individual abundances can be provided by first giving how many of them are changed and then
# for each of them their atomic number followed by the absolute abundance on the same line
#
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
26  $Feabu
###########
# if xifix true, fixed microturbulence is read from next line (km/s)
# otherwise the value(s) are read from the model atmosphere.
#
'XIFIX:' 'F'
$TURBVEL
EOF

########################################################################
set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}-NLTE-windows.spec
set result     = ${MODEL}${SUFFIX}

time ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/exec/bsyn_lu <<EOF
###########
# Use NLTE if true. Source function is computed with departure coefficients 
# from departure coefficient file for the atom in model atom file, if they 
# are provided. 
#
'NLTE :'          '.true.'
###########
# file containing NLTE information (species and associated files)
#
'NLTEINFOFILE:'  'DATA/SPECIES_LTE_NLTE.dat'
#
###########
# if present these files will be used to compute the spectrum in a number of 
# windows given by SEGMENTSFILE. 
# Comment out if not needed
#
'SEGMENTSFILE:'     '${dpath}/uves_giant_Fe-seg.txt'
#
# spectral resolution  to be used in these windows. 
# If not specified, a default value of 500000 is used
#
'RESOLUTION:'     '300000.'
###########
# spectral interval in the case of a single wavelength interval, i.e. no 
# segmentsfile.
# min and max lambda for the calculations and constant wavelength step
#
'LAMBDA_MIN:'     '${lam_min}'
'LAMBDA_MAX:'     '${lam_max}'
'LAMBDA_STEP:'    '${deltalam}'
###########
# Intensity / Flux 
#
'INTENSITY/FLUX:' 'Flux'
###########
# for eqwidt only: if true, iterate abundance until observed eqwivalent 
# widths are matched. NOT IMPLEMENTED IN THIS VERSION
#
'ABFIND        :' '.false.'
###########
# file containing continuous opacity at all model points. Can be computed with babsma.f
# or be a .opa file from the MARCS web site.
#
'MODELOPAC:' 'contopac/${MODEL}opac'
###########
# output file containing spectrum, or equivalent widths
#
'RESULTFILE :' 'syntspec/${result}'
###########
# chemical composition
#
# in bsyn.f, in addition to what is explained above for babsma.f,
# isotopic fractions can be set, e.g. 
# 12.012 0.9
# 12.013 0.1
# to set 90% of 12C and 10% of 13C.
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
26  $Feabu
'ISOTOPES : ' '2'
12.012 0.9
12.013 0.1
###########
# line lists. First how many there are, and then the list of lists
#
'NFILES   :' '2'
TEST-data/nlte_linelist_test.txt
DATA/Hlinedata
###########
# spherical (T) or plane-parallel (F) radiative transfer. 
# If spherical, a few more parameters are read
# DO NOT CHANGE THESE PARAMETERS UNLESS YOU REALLY KNOW WHAT YOU ARE DOING.
#
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF
