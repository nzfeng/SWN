set(GUROBI_CUSTOM /Library/gurobi1000/macos_universal2) # replace with path to Gurobi on your computer

find_path(GUROBI_INCLUDE_DIRS
    NAMES gurobi_c++.h
    HINTS ${GUROBI_CUSTOM}/include ${GUROBI_HOME} ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES include)

find_library(GUROBI_LIBRARY
    NAMES gurobi gurobi100
    HINTS ${GUROBI_CUSTOM}/lib ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib)


if(MSVC)
    set(MSVC_YEAR "2017")
    
    if(MT)
        set(M_FLAG "mt")
    else()
        set(M_FLAG "md")
    endif()
    
    find_library(GUROBI_CXX_LIBRARY
        NAMES gurobi_c++${M_FLAG}${MSVC_YEAR}
        HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
        PATH_SUFFIXES lib)
    find_library(GUROBI_CXX_DEBUG_LIBRARY
        NAMES gurobi_c++${M_FLAG}d${MSVC_YEAR}
        HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
        PATH_SUFFIXES lib)
else()
    find_library(GUROBI_CXX_LIBRARY
        NAMES gurobi_c++ libgurobi_c++.a
        HINTS ${GUROBI_CUSTOM}/lib ${GUROBI_DIR} ${GUROBI_HOME} $ENV{GUROBI_HOME}
        PATH_SUFFIXES lib)
    message(STATUS "GUROBI_CXX_LIBRARY: ${GUROBI_CXX_LIBRARY}")
    set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_CXX_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY)
