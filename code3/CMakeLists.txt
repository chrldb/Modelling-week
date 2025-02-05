cmake_minimum_required(VERSION 3.10)
project(ElectricFieldSim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ✅ Trouver OpenMP (pour parallélisation)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    add_compile_options(${OpenMP_CXX_FLAGS})
    add_link_options(${OpenMP_CXX_FLAGS})
endif()

# ✅ Compilation optimisée (-O3 pour rapidité)
add_executable(ElectricFieldSim main.cpp)

# ✅ Linker avec OpenMP (si utilisé)
target_link_libraries(ElectricFieldSim OpenMP::OpenMP_CXX)

# ✅ Supprimer les warnings inutiles
add_compile_options(-Wno-deprecated-declarations)