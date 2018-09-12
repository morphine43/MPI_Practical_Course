cmake_minimum_required( VERSION 2.8 )

message( STATUS "MPI_Practical_Course" )

find_package( MPI REQUIRED )
if( MPI_FOUND )
    include_directories( ${MPI_INCLUDE_PATH} )
endif( MPI_FOUND )

############ test functionality ############
if ( NOT MSVC )
    find_package( OpenCV REQUIRED )
    if( OpenCV_FOUND )
        message(STATUS "OpenCV - Found")
        include_directories( ${OpenCV_INCLUDE_DIRS} )
    endif( OpenCV_FOUND )
endif( NOT MSVC )
############################################

add_subdirectory(modules)
