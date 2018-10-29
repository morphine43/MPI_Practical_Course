# Install script for directory: C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files/Project")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/Perov_Dima_task1_SumElVect/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/vikhrev_array_sum/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/Yakovlev_Pavel_mul_vect/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/Vdovin_Eugene_task1_NumWords/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/Repin_Vladimir_Wrong_Order_Count/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/Gaidaichuk_Yuri_task1_Arrays_Matching/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/Druzhinin_Alexei_task1_MinElVect/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/Zolotareva_Olesya_task1_MiddleValueVect/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/ndanilov_rows_matrix_sum/cmake_install.cmake")
  include("C:/Users/Agricola/Documents/UNN/MPI/task2_1/MPI_Practical_Course/modules/task_1/Yunin_Ivan_maximum_elem_vector/cmake_install.cmake")

endif()

