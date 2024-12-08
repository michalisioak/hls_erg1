cmake_minimum_required(VERSION 3.14)


project(erg1 VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(AC_MATH_DIR "./ac_math") 
set(AC_TYPES_DIR "./ac_types") 

set(SOURCE_FILES src/main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_compile_definitions(${PROJECT_NAME} PRIVATE TESTING)

target_include_directories(${PROJECT_NAME} PRIVATE ${AC_MATH_DIR}/include)
target_include_directories(${PROJECT_NAME} PRIVATE ${AC_TYPES_DIR}/include)

# target_link_libraries(${PROJECT_NAME} ac_math)
# target_link_libraries(${PROJECT_NAME} ac_math)


target_compile_options(${PROJECT_NAME} PRIVATE -Wall -O0 -Wextra -Wno-deprecated-copy -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-parameter
)
