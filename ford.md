project: sqpopt
src_dir: ./src
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/sqpopt
summary: Modern Fortran SQP OPTimizer
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
html_template_dir: ./ford/templates
css: ./ford/user.css
source: true
graph: true
extra_mods: fmin_module:https://github.com/jacobwilliams/fmin
            iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            lusol:https://github.com/jacobwilliams/lusol
            lusol_precision:https://github.com/jacobwilliams/lusol
            lsqr_module:https://github.com/jacobwilliams/LSQR

{!README.md!}
