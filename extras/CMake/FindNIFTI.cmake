# ZLIB is required
# ================
FIND_PACKAGE( ZLIB REQUIRED )


# Search for nifticlib headers and libraries
# ==========================================
FIND_PATH( NIFTI_INCLUDE_DIR
	nifti1.h
	/usr/local/include/nifti /usr/include/nifti
)

FIND_LIBRARY( NIFTI_BASE_LIBRARY
	NAMES libniftiio${CMAKE_SHARED_LIBRARY_SUFFIX} libniftiio${CMAKE_STATIC_LIBRARY_SUFFIX}
)

FIND_LIBRARY( NIFTI_ZNZ_LIBRARY
	NAMES libznz${CMAKE_SHARED_LIBRARY_SUFFIX} libznz${CMAKE_STATIC_LIBRARY_SUFFIX}
)


# If everything is ok, set variable for the package
# =================================================
IF( NIFTI_INCLUDE_DIR AND NIFTI_BASE_LIBRARY AND NIFTI_ZNZ_LIBRARY )

   SET( NIFTI_FOUND TRUE )
   SET( NIFTI_INCLUDE_DIRS
      ${NIFTI_INCLUDE_DIR}
      ${ZLIB_INCLUDE_DIRS}
   )
   SET( NIFTI_LIBRARIES
      ${NIFTI_BASE_LIBRARY}
      ${NIFTI_ZNZ_LIBRARY}
      ${ZLIB_LIBRARIES}
   )
   IF( NOT NIFTI_FIND_QUIETLY )
      MESSAGE(STATUS "Found NIFTI: ${NIFTI_LIBRARIES}")
   ENDIF (NOT NIFTI_FIND_QUIETLY)

ELSE()

   SET( NIFTI_FOUND FALSE )
   IF( NIFTI_FIND_REQUIRED )
      MESSAGE(FATAL_ERROR "Could not find NIFTI (which is required)")
   ENDIF()

ENDIF()
