#!/bin/csh -f

#
# This script uses a line list with measured equivalent widths (and error bars)
# to derive the line by line abundance of the element. The output is a list of
# lines with an eqw measurement and the corresponding abundance (bracketed 
# by the abundance error.
# each line in the list is treated as a single line. Blends or hyperfine
# structure cannot be treated properly.
# As the molecular equilibrium is only calculated once, for the input abundances
# specified in the script, it is highly recommended to rerun the script with
# modified abundances close to the abundances deduced from the eqw.
#

date
set mpath=models

foreach MODEL (5777g4.44z+0.00a+0.00t01-ref107857)

#set Cabu = 8.656
#set Nabu = 7.78
#set Oabu = 8.66
 
set lam_min    = '6700'
set lam_max    = '6720'

set deltalam   = '0.1'
set METALLIC   = '     0.000'
set TURBVEL    = '1.0'
set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}.findabu
set result     = ${MODEL}${SUFFIX}

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

../exec-v19.2/eqwidt_lu <<EOF
'LAMBDA_MIN:'     '${lam_min}'
'LAMBDA_MAX:'     '${lam_max}'
'LAMBDA_STEP:'    '${deltalam}'
'INTENSITY/FLUX:' 'Flux'
'COS(THETA)    :' '1.00'
'ABFIND        :' '.true.'
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
'NFILES   :' '1'
linelists/vald_eqw.list
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
EOF
########################################################################
date
end
