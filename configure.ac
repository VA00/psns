#    -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([libpsns],[0.12],[andrzej.odrzywolek@uj.edu.pl],[],[2024-11-12])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror gnu])
AC_CONFIG_MACRO_DIR([m4])
        
AC_CONFIG_FILES([                Makefile 
                             src/Makefile
])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AM_PROG_AR
LT_INIT

# Checks for libraries.




AC_CHECK_LIB([m], [pow, sqrt, log10, exp])
AC_CHECK_LIB([cuba], [Cuhre, Divonne, Vegas, Suave],[],
[
echo "Cuba library not found!"
exit -1
]
)
AC_CHECK_LIB([nse],[NSE_enum],[],
[
echo "libnse not found!"
exit -1
]
)
AC_CHECK_LIB([interp],[trilinear_interp],[],
[
echo "libinterp not found!"
exit -1
])

AC_CHECK_LIB([fermidirac], [Ffermi], [], [
    echo "libfermidirac not found!"
    exit -1
], [-lflint])

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

# Checks for header files.

AC_CHECK_HEADER([cuba.h], [AC_DEFINE([HAVE_CUBA_H])], [AC_MSG_ERROR([Cuba multidimensional integration library missing:  Computing exact pair annihilation spectrum impossible! Download and install from www.feynarts.de/cuba/ and/or add information to ./configure from command line.])])

AC_CHECK_HEADER([nse.h], [AC_DEFINE([HAVE_NSE_H])], [AC_MSG_ERROR([NSE lib header not found!])])

AC_CHECK_HEADER([interp.h], [AC_DEFINE([HAVE_INTERP_H])], [AC_MSG_ERROR([Trilinear interpolation header interp.h not found!])])

AC_CHECK_HEADER([fermidirac.h],     [AC_DEFINE([HAVE_FERMIDIRAC_H])],     [AC_MSG_ERROR([Fermi-Dirac functions header fermidirac.h not found!      ])])


AC_CHECK_HEADERS([stdio.h stdlib.h string.h unistd.h math.h nse.h interp.h cuba.h fermidirac.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL


AC_OUTPUT
