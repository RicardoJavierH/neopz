function(enable_lapack target)
    if(USING_MKL)
        set(BLA_VENDOR "Intel10_64lp") #tries to find MKL version
    endif()
    find_package(LAPACK REQUIRED)
    target_link_libraries(${target} PRIVATE ${LAPACK_LIBRARIES} )
    target_compile_definitions(${target} PRIVATE USING_LAPACK)
    target_compile_definitions(${target} PRIVATE USING_BLAS)
    foreach(LA_LIB LAPACK_LIBRARIES)
        if(${LA_LIB} MATCHES ".*mkl.*")
            message("Found LAPACK from MKL")
            target_compile_definitions(${target} PRIVATE MKLLAPACK)
            target_compile_definitions(${target} PRIVATE MKLBLAS)
            get_filename_component(LA_LIB_DIR ${LA_LIB} DIRECTORY)
            include(cmake/FindIncludeFromMKL.cmake)
            find_include_from_mkl(LAPACK_INCLUDE_DIR ${BL_LIB_DIR} mkl_lapacke.h)
            if(NOT LAPACK_INCLUDE_DIR)
              message(FATAL_ERROR "Could not find LAPACK include directory")
            else()
              message(STATUS "LAPACK include dirs: ${LAPACK_INCLUDE_DIR}")
              target_include_directories(${target} PRIVATE ${LAPACK_INCLUDE_DIR})
              mark_as_advanced(LAPACK_INCLUDE_DIR)
            endif()
            break()
        endif()
    endforeach()
endfunction()