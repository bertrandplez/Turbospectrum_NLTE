#!/bin/csh -f

date
set mpath = ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/COM/TEST-data/
set dpath = ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/COM/TEST-data/

set MODEL = atmos.sun_marcs_t5777_4.44_0.00_vmic1_new

set Caabu = 6.34

#set Cabu = 8.656
#set Nabu = 7.78
set Oabu = 8.66
 
set lam_min    = '4800.'
set lam_max    = '6800.'

set deltalam   = '0.006'
set METALLIC   = '     0.000'
set TURBVEL = 1.0

#
# ABUNDANCES FROM THE MODEL ARE NOT USED !!!

time ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/exec/babsma_lu <<EOF
###########
# into the absorption part. This happens in babsma.f only. SHOULD BE IMPLEMENTED IN BSYN ALSO!
#
'PURE-LTE  :'  '.false.'
###########
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
# in bsyn.f isotopic fractions can be set, e.g. 
# 12.012 0.9
# 12.013 0.1
# to set 90% of 12C and 10% of 13C.
#
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
20  $Caabu
###########
# if xifix true, fixed microturbulence is read from next line (km/s)
# otherwise the value(s) are read from the model atmosphere.
#
'XIFIX:' 'F'
$TURBVEL
EOF

########################################################################
set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}-LTE-windows.spec
set result     = ${MODEL}${SUFFIX}

time ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/exec/bsyn_lu <<EOF
###########
# if PURE-LTE is true, the scattering part of the continuum opacity is added
# into the absorption part. This happens in babsma.f only. SHOULD BE IMPLEMENTED IN BSYN ALSO!
# Nothing happens if it is set to true in bsyn.f only.
#
'PURE-LTE  :'  '.false.'
###########
# Use NLTE if true. Source function is computed with departure coefficients 
# from departure coefficient file for the atom in model atom file, if they 
# are provided. Line list must include the same species with a 'NLTE' flag.
#
'NLTE :'          '.false.'
###########
# file containing model atom (levels and energies)
#
'MODELATOMFILE:'  ' '
###########
# departure coefficients file. It must have been computed for the same 
# model atomsphere used here.
#
'DEPARTUREFILE:'  ' '
###########
# if present these files will be used to compute the spectrum in a number of 
# windows givent by LINEMASKFILE and CONTMASKFILE, using lines and continuum
# opacities collected in windows given by SEGMENTSFILE. 
# Comment out if not needed
#
'CONTMASKFILE:'     '${dpath}/uves_giant_Fe-cmask.txt'
'LINEMASKFILE:'     '${dpath}/uves_giant_Fe-lmask.txt'
'SEGMENTSFILE:'     '${dpath}/uves_giant_Fe-seg.txt'
###########
# spectral interval
#
'LAMBDA_MIN:'     '${lam_min}'
'LAMBDA_MAX:'     '${lam_max}'
'LAMBDA_STEP:'    '${deltalam}'
###########
# Intensity / Flux followed by angle at which intensity should be computed.
#
'INTENSITY/FLUX:' 'Flux'
'COS(THETA)    :' '1.00'
###########
# for eqwidt only: if true, iterate abundance until observed eqwivalent 
# widths are matched. 
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
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
20  $Caabu
'ISOTOPES : ' '0'
###########
# line lists. First how many there are, and then the list
#
'NFILES   :' '2'
DATA/Hlinedata
~/Documents/GitHub/Turbospectrum/Turbospectrum2020/COM/TEST-data/lte_linelist_test.txt
#'NFILES   :' '1'
#~/Documents/GitHub/Turbospectrum/Turbospectrum2020/COM/TEST-data/CaI_lte_linelist.dat
###########
# spherical or plane-parallel radiative transfer. If spherical, a few more parameters are read
#
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF
########################################################################
set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}-LTE.spec
set result     = ${MODEL}${SUFFIX}

time ~/Documents/GitHub/Turbospectrum/Turbospectrum2020/exec/bsyn_lu <<EOF
###########
# if PURE-LTE is true, the scattering part of the continuum opacity is added
# into the absorption part. This happens in babsma.f only. SHOULD BE IMPLEMENTED IN BSYN ALSO!
# Nothing happens if it is set to true in bsyn.f only.
#
'PURE-LTE  :'  '.false.'
###########
# Use NLTE if true. Source function is computed with departure coefficients 
# from departure coefficient file for the atom in model atom file, if they 
# are provided. Line list must include the same species with a 'NLTE' flag.
#
'NLTE :'          '.false.'
###########
# file containing model atom (levels and energies)
#
'MODELATOMFILE:'  ' '
###########
# departure coefficients file. It must have been computed for the same 
# model atomsphere used here.
#
'DEPARTUREFILE:'  ' '
###########
# if present these files will be used to compute the spectrum in a number of 
# windows givent by LINEMASKFILE and CONTMASKFILE, using lines and continuum
# opacities collected in windows given by SEGMENTSFILE. 
# Comment out if not needed
#
#'CONTMASKFILE:'     '${dpath}/uves_giant_Fe-cmask.txt'
#'LINEMASKFILE:'     '${dpath}/uves_giant_Fe-lmask.txt'
#'SEGMENTSFILE:'     '${dpath}/uves_giant_Fe-seg.txt'
###########
# spectral interval
#
'LAMBDA_MIN:'     '${lam_min}'
'LAMBDA_MAX:'     '${lam_max}'
'LAMBDA_STEP:'    '${deltalam}'
###########
# Intensity / Flux followed by angle at which intensity should be computed.
#
'INTENSITY/FLUX:' 'Flux'
'COS(THETA)    :' '1.00'
###########
# for eqwidt only: if true, iterate abundance until observed eqwivalent 
# widths are matched. 
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
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
20  $Caabu
'ISOTOPES : ' '0'
###########
# line lists. First how many there are, and then the list
#
'NFILES   :' '2'
DATA/Hlinedata
~/Documents/GitHub/Turbospectrum/Turbospectrum2020/COM/TEST-data/lte_linelist_test.txt
#'NFILES   :' '1'
#~/Documents/GitHub/Turbospectrum/Turbospectrum2020/COM/TEST-data/CaI_lte_linelist.dat
###########
# spherical or plane-parallel radiative transfer. If spherical, a few more parameters are read
#
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF
