# Matlab binding to CPUfit

cmake_minimum_required(VERSION 3.7)

find_package( MATLAB REQUIRED MX_LIBRARY REQUIRED)

if(NOT Matlab_MX_LIBRARY)
	return()
endif()

set( TestHeaders
)

set( TestSources
	cminpackMex.cpp
)

	
add_library( cminpackMex SHARED
	${TestHeaders}
	${TestSources}
)

target_link_libraries( cminpackMex cminpack ${Matlab_LIBRARIES} )
set_property( TARGET cminpackMex
	PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}" )

include_directories( ${Matlab_INCLUDE_DIRS} )

if( WIN32 )
	if( CMAKE_CL_64 )
		set( suffix .mexw64 )
	else()
		set( suffix .mexw32 )
	endif()
	set_target_properties( cminpackMex
		PROPERTIES
			OUTPUT_NAME cminpackMex
			SUFFIX ${suffix}
	)
	SET(CMAKE_SHARED_LINKER_FLAGS "/export:mexFunction")
endif()
