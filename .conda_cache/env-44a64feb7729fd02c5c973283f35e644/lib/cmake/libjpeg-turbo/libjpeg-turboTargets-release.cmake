#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libjpeg-turbo::jpeg" for configuration "Release"
set_property(TARGET libjpeg-turbo::jpeg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::jpeg PROPERTIES
  IMPORTED_LOCATION_RELEASE "/Users/saba/Documents/GitHub/PangenomeAnalysis/.conda_cache/env-44a64feb7729fd02c5c973283f35e644/lib/libjpeg.8.3.2.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libjpeg.8.dylib"
  )

list(APPEND _cmake_import_check_targets libjpeg-turbo::jpeg )
list(APPEND _cmake_import_check_files_for_libjpeg-turbo::jpeg "/Users/saba/Documents/GitHub/PangenomeAnalysis/.conda_cache/env-44a64feb7729fd02c5c973283f35e644/lib/libjpeg.8.3.2.dylib" )

# Import target "libjpeg-turbo::turbojpeg" for configuration "Release"
set_property(TARGET libjpeg-turbo::turbojpeg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::turbojpeg PROPERTIES
  IMPORTED_LOCATION_RELEASE "/Users/saba/Documents/GitHub/PangenomeAnalysis/.conda_cache/env-44a64feb7729fd02c5c973283f35e644/lib/libturbojpeg.0.5.0.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libturbojpeg.0.dylib"
  )

list(APPEND _cmake_import_check_targets libjpeg-turbo::turbojpeg )
list(APPEND _cmake_import_check_files_for_libjpeg-turbo::turbojpeg "/Users/saba/Documents/GitHub/PangenomeAnalysis/.conda_cache/env-44a64feb7729fd02c5c973283f35e644/lib/libturbojpeg.0.5.0.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
