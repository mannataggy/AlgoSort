#cmake_minimum_required(VERSION 3.22)
#project(sorting-visualizer)

#set(SFML_DIR "C:\SFML-2.5.1\lib\cmake\SFML\SFMLConfig.cmake")

#set(SFML_DIR "file:///C:/SFML-2.5.1/lib/cmake")

#find_package(SFML 2.5 COMPONENTS graphics REQUIRED)

#find_package(SFML COMPONENTS graphics window REQUIRED)

#add_executable(sorting-visualizer src/main.cpp src/Sortable.cpp src/SortAlgorithms.cpp src/SortController.cpp src/Utils.cpp)

#target_link_libraries(sorting-visualizer sfml-graphics sfml-window)

cmake_minimum_required(VERSION 3.5)

# Set the policy to use <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

project(sorting-visualizer)

set(SFML_DIR "C:\\SFML-2.5.1\\lib\\cmake\\SFML")

# Set the path to your SFML installation
set(SFML_ROOT "C:/SFML-2.5.1")

list(APPEND CMAKE_MODULE_PATH "C:\\SFML-2.5.1\\lib\\cmake\\SFML")

#find_package(SFML 2.5 COMPONENTS graphics REQUIRED)
find_package(SFML COMPONENTS graphics window REQUIRED)

# Add your source files
add_executable(sorting-visualizer src/main.cpp src/Sortable.cpp src/SortAlgorithms.cpp src/SortController.cpp src/Utils.cpp)

include_directories("C:/SFML/include")
link_directories("C:/SFML/lib")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/Modules")

# Link your executable with the SFML libraries
target_link_libraries(sorting-visualizer sfml-graphics)

# Copy the required SFML DLLs to the build directory
add_custom_command(TARGET sorting-visualizer POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${SFML_ROOT}/bin" $<TARGET_FILE_DIR:sorting-visualizer>)
