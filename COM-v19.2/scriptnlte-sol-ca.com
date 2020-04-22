#!/bin/csh -f

date
set mpath = TEST-data/

set MODEL = atmos.sun_marcs_t5777_4.44_0.00_vmic1_new

#set Cabu = 8.656
#set Nabu = 7.78
set Oabu = 8.66
 
set lam_min    = '5330.'
set lam_max    = '5630.'

set deltalam   = '0.01'
set METALLIC   = '     0.000'
set TURBVEL    = '1.0'


../exec-v19.2/babsma_lu << EOF
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
'INDIVIDUAL ABUNDANCES:'   '0'
'XIFIX:' 'T'
$TURBVEL
EOF

########################################################################
set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}-v19.2-NLTE.spec
set result     = ${MODEL}${SUFFIX}

../exec-v19.2/bsyn_lu <<EOF
'PURE-LTE  :'  '.false.'
'NLTE :'          '.true.'
'MODELATOMFILE:'  'TEST-data/atom.caNew'
'DEPARTUREFILE:'  'TEST-data/caNew_output_4TS.csv'
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
'INDIVIDUAL ABUNDANCES:'   '0'
'ISOTOPES : ' '0'
'NFILES   :' '1'
TEST-data/CaI_nlte_linelist.dat
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF

########################################################################
set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}-v19.2-LTE.spec
set result     = ${MODEL}${SUFFIX}

../exec-v19.2/bsyn_lu <<EOF
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
'INDIVIDUAL ABUNDANCES:'   '0'
'ISOTOPES : ' '0'
'NFILES   :' '1'
TEST-data/CaI_lte_linelist.dat
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF
########################################################################
