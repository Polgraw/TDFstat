#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "sleef::sleef" for configuration "Release"
set_property(TARGET sleef::sleef APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(sleef::sleef PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libsleef.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS sleef::sleef )
list(APPEND _IMPORT_CHECK_FILES_FOR_sleef::sleef "${_IMPORT_PREFIX}/lib64/libsleef.a" )

# Import target "sleef::sleefgnuabi" for configuration "Release"
set_property(TARGET sleef::sleefgnuabi APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(sleef::sleefgnuabi PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libsleefgnuabi.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS sleef::sleefgnuabi )
list(APPEND _IMPORT_CHECK_FILES_FOR_sleef::sleefgnuabi "${_IMPORT_PREFIX}/lib64/libsleefgnuabi.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
