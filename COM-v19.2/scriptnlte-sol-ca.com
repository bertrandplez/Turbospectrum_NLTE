#!/bin/csh -f

date
set mpath = ~/Documents/GitHub/Turbospectrum/Turbospectrum2019.2/COM-v19.2/TEST-data/

set MODEL = atmos.sun_marcs_t5777_4.44_0.00_vmic1_new

set Caabu = 6.34

#set Cabu = 8.656
#set Nabu = 7.78
set Oabu = 8.66
 
set lam_min    = '5349.'
set lam_max    = '5350.'

set deltalam   = '0.005'
set METALLIC   = '     0.000'

#
# ABUNDANCES FROM THE MODEL ARE NOT USED !!!

~/Documents/GitHub/Turbospectrum/Turbospectrum2019.2/exec-v19.2/babsma_lu << EOF
'PURE-LTE  :'  '.false.'
'LAMBDA_MIN:'  '${lam_min}'
'LAMBDA_MAX:'  '${lam_max}'
'LAMBDA_STEP:' '${deltalam}'
'MODELINPUT:' '$mpath/${MODEL}'
'MARCS-FILE:' '.false.'
'MODELOPAC:' 'contopac/${MODEL}opac'
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
20  $Caabu
'XIFIX:' 'F'
$TURBVEL
EOF

########################################################################
set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}-v19.2-LTE.spec
set result     = ${MODEL}${SUFFIX}

~/Documents/GitHub/Turbospectrum/Turbospectrum2019.2/exec-v19.2/bsyn_lu <<EOF
'PURE-LTE  :'  '.false.'
'NLTE :'          '.false.'
'MODELATOMFILE:'  ' '
'DEPARTUREFILE:'  ' '
'LAMBDA_MIN:'     '${lam_min}'
'LAMBDA_MAX:'     '${lam_max}'
'LAMBDA_STEP:'    '${deltalam}'
'INTENSITY/FLUX:' 'Flux'
'COS(THETA)    :' '1.00'
'ABFIND        :' '.false.'
'MODELOPAC:' 'contopac/${MODEL}opac'
'RESULTFILE :' 'syntspec/${result}'
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
20  $Caabu
'ISOTOPES : ' '0'
'NFILES   :' '1'
~/Documents/GitHub/Turbospectrum/Turbospectrum2019.2/COM-v19.2/TEST-data/CaI_lte_linelist.dat
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF
########################################################################
set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}-v19.2-NLTE.spec
set result     = ${MODEL}${SUFFIX}

~/Documents/GitHub/Turbospectrum/Turbospectrum2019.2/exec-v19.2/bsyn_lu <<EOF
'PURE-LTE  :'  '.false.'
'NLTE :'          '.true.'
'MODELATOMFILE:'  '~/Documents/GitHub/Turbospectrum/Turbospectrum2019.2/COM-v19.2/TEST-data/atom.caNew'
'DEPARTUREFILE:'  '~/Documents/GitHub/Turbospectrum/Turbospectrum2019.2/COM-v19.2/TEST-data/sun_marcs_t5777_4.44_0.00_vmic1_new_A_Ca_6.34.bin'
'LAMBDA_MIN:'     '${lam_min}'
'LAMBDA_MAX:'     '${lam_max}'
'LAMBDA_STEP:'    '${deltalam}'
'INTENSITY/FLUX:' 'Flux'
'COS(THETA)    :' '1.00'
'ABFIND        :' '.false.'
'MODELOPAC:' 'contopac/${MODEL}opac'
'RESULTFILE :' 'syntspec/${result}'
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
20  $Caabu
'ISOTOPES : ' '0'
'NFILES   :' '1'
~/Documents/GitHub/Turbospectrum/Turbospectrum2019.2/COM-v19.2/TEST-data/CaI_nlte_linelist.dat
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF

########################################################################
