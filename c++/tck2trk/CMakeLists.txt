PROJECT( tck2trk )
ADD_EXECUTABLE( ${PROJECT_NAME} main.cxx )


########################  FIND PACKAGES USED  #######################################
FIND_PACKAGE( BLITZ REQUIRED )
FIND_PACKAGE( NIFTI REQUIRED )


############################  INCLUDE / LINK  #######################################
INCLUDE_DIRECTORIES( ${NIFTI_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${BLITZ_INCLUDE_DIR} )

TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
	${BLITZ_LIBRARIES}
	${NIFTI_LIBRARIES}
)

###################################  INSTALL  #######################################
INSTALL( TARGETS ${PROJECT_NAME}
	RUNTIME DESTINATION "bin"
)

