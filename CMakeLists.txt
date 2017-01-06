# The name of our project is "CMINPACK". CMakeLists files in this project can
# refer to the root source directory of the project as ${CMINPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${CMINPACK_BINARY_DIR}.



#include(${PROJECT_SOURCE_DIR}/cminpack-1.3.4/cmake/cminpack_utils.cmake)
# Set version and OS-specific settings
#set(CMINPACK_VERSION 1.3.4 CACHE STRING "CMinpack version")
#set(CMINPACK_SOVERSION 1 CACHE STRING "CMinpack API version")
#DISSECT_VERSION()


#set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../build)

add_subdirectory (cmake)

project( cminpack )

set (cminpack_srcs
  cminpack.h cminpackP.h
  chkder.c  enorm.c   hybrd1.c  hybrj.c   lmdif1.c  lmstr1.c  qrfac.c   r1updt.c
  dogleg.c  fdjac1.c  hybrd.c   lmder1.c  lmdif.c   lmstr.c   qrsolv.c  rwupdt.c
  dpmpar.c  fdjac2.c  hybrj1.c  lmder.c   lmpar.c   qform.c   r1mpyq.c  covar.c covar1.c
  minpack.h
  chkder_.c enorm_.c  hybrd1_.c hybrj_.c  lmdif1_.c lmstr1_.c qrfac_.c  r1updt_.c
  dogleg_.c fdjac1_.c hybrd_.c  lmder1_.c lmdif_.c  lmstr_.c  qrsolv_.c rwupdt_.c
  dpmpar_.c fdjac2_.c hybrj1_.c lmder_.c  lmpar_.c  qform_.c  r1mpyq_.c covar_.c
  )
set (cminpack_hdrs
    cminpack.h minpack.h)

add_library (cminpack SHARED ${cminpack_srcs})

set_property( TARGET cminpack
	PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}" )
	
#set_target_properties(cminpack PROPERTIES VERSION ${CMINPACK_VERSION} SOVERSION ${CMINPACK_SOVERSION})

add_subdirectory (matlab)