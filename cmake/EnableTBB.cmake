function(enable_tbb target)
    find_package(TBB REQUIRED)
    if (NOT TBB_FOUND)
        message(FATAL_ERROR "Could not find TBB. If TBB is not needed, "
                "configure the project with -DUSING_TBB=OFF")
    endif()
    target_link_libraries(${target} PRIVATE TBB::tbb)
    target_include_directories(${target} PRIVATE ${TBB_INCLUDE_DIRS})
endfunction()
