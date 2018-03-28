# Try to find the libssh library
#  LIBSSH_FOUND - system has libssh lib
#  LIBSSH_INCLUDE_DIRS - the include directories needed
#  LIBSSH_LIBRARIES - libraries needed

FIND_PATH(LIBSSH_INCLUDE_DIR NAMES libssh/libssh.h)
FIND_LIBRARY(LIBSSH_LIBRARY NAMES ssh libssh)
SET(LIBSSH_INCLUDE_DIRS ${LIBSSH_INCLUDE_DIR})
SET(LIBSSH_LIBRARIES ${LIBSSH_LIBRARY})
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBSSH DEFAULT_MSG LIBSSH_LIBRARY LIBSSH_INCLUDE_DIR)
MARK_AS_ADVANCED(LIBSSH_INCLUDE_DIR LIBSSH_LIBRARY)