#!/bin/csh -f
# compute intensity spectrum for the Sun, at 2 different angles


date
set mpath=models

foreach MODEL (5777g4.44z+0.00a+0.00t01-ref107857)

#set Cabu = 8.656
#set Nabu = 7.78
#set Oabu = 8.66
 
set lam_min    = '6700'
set lam_max    = '6720'

set deltalam   = '0.01'
set METALLIC   = '     0.000'
set TURBVEL    = '1.0'

#
# ABUNDANCES FROM THE MODEL ARE NOT USED !!!

../exec-v19.2/babsma_lu << EOF
'LAMBDA_MIN:'  '${lam_min}'
'LAMBDA_MAX:'  '${lam_max}'
'LAMBDA_STEP:' '${deltalam}'
'MODELINPUT:' '$mpath/${MODEL}'
'MARCS-FILE:' '.true.'
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

foreach MU (1.0 0.3)
# mu=cos(theta); mu=1.0 is disk center

set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}_mu${MU}.intensity
set result     = ${MODEL}${SUFFIX}

../exec-v19.2/bsyn_lu <<EOF
'LAMBDA_MIN:'     '${lam_min}'
'LAMBDA_MAX:'     '${lam_max}'
'LAMBDA_STEP:'    '${deltalam}'
'INTENSITY/FLUX:' 'Intensity'
'COS(THETA)    :' '${MU}'
'ABFIND        :' '.false.'
'MODELOPAC:' 'contopac/${MODEL}opac'
'RESULTFILE :' 'syntspec/${result}'
'METALLICITY:'    '${METALLIC}'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
3  1.05
'ISOTOPES : ' '2'
3.006  0.075
3.007  0.925
'NFILES   :' '2'
DATA/Hlinedata
linelists/vald-6700-6720.list
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF

end

########################################################################
date
end
