AC_DEFUN([AC_PROG_DOXYGEN],
[
AC_ARG_ENABLE([doxygen], [  --enable-doxygen        enable documentation generation with doxygen (auto)])                          

if test "x$enable_doxygen" = xno; then
 enable_doc=no	
else 
 AC_CHECK_PROG([DOXYGEN],[doxygen],[doxygen])
 if test x$DOXYGEN = x; then
  if test "x$enable_doxygen" = xyes; then
   AC_MSG_ERROR([doxygen not found])
  fi
  enable_doc=no
 else
  enable_doc=yes
  AC_CHECK_PROG([PDFLATEX],[pdflatex],[pdflatex])
 fi
fi

AM_CONDITIONAL([DOC],[test x$enable_doc = xyes])

if tes x$PDFLATEX = x; then
 AC_MSG_ERROR([pdflatex not found])
fi

])

