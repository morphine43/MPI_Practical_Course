cmake_minimum_required( VERSION 2.8 )

message( STATUS "MPI_Practical_Course" )

find_package( MPI REQUIRED )
if( MPI_FOUND )
    include_directories( ${MPI_INCLUDE_PATH} )
endif( MPI_FOUND )

add_subdirectory(modules)
