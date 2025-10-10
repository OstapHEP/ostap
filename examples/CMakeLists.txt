           
install ( DIRECTORY      ${CMAKE_CURRENT_SOURCE_DIR}  
          DESTINATION    .
          COMPONENT      examples 
          FILES_MATCHING PATTERN "*.py"
          PERMISSIONS    OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install ( DIRECTORY      ${CMAKE_CURRENT_SOURCE_DIR}  
          DESTINATION    .  
          COMPONENT      examples 
          FILES_MATCHING PATTERN "*.C"
          PERMISSIONS    OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install ( DIRECTORY      ${CMAKE_CURRENT_SOURCE_DIR}  
          DESTINATION    .  
          COMPONENT      examples 
          FILES_MATCHING PATTERN "*.cpp"
          PERMISSIONS    OWNER_READ GROUP_READ  WORLD_READ )
   
install ( DIRECTORY      ${CMAKE_CURRENT_SOURCE_DIR}  
          DESTINATION    .
          COMPONENT      examples 
          FILES_MATCHING PATTERN "*.cxx"
          PERMISSIONS    OWNER_READ GROUP_READ  WORLD_READ )
   
enable_testing() 

file(GLOB pyfiles "${CMAKE_INSTALL_PREFIX}/examples/*/*.py")
foreach(filename ${pyfiles})
   get_filename_component(fname ${filename} NAME_WE)
   string(CONCAT fname ostap.examples_py_ ${fname})
   add_test (NAME ${fname} COMMAND python     ${filename} )
   set_tests_properties(${fname} PROPERTIES TIMEOUT 5000 ) 
endforeach() 

file(GLOB cfiles "${CMAKE_INSTALL_PREFIX}/examples/*/*.C")
foreach(cfilename ${cfiles})
   get_filename_component(cfname ${cfilename} NAME_WE)
   string(CONCAT cfname ostap.examples_C_ ${cfname})
   add_test (NAME ${cfname}C WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary2COMMAND root.exe ${cfilename} -b -q )
   set_tests_properties(${cfname}C PROPERTIES TIMEOUT 5000 )
endforeach() 

file(GLOB cppfiles "${CMAKE_INSTALL_PREFIX}/examples/*/*.cpp")
foreach(cppfilename ${cppfiles})
   get_filename_component(cppfname ${cppfilename} NAME_WE)
   string(CONCAT cppfname ostap.examples_cpp_ ${cppfname})
   add_test (NAME ${cppfname}CPP WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary2 COMMAND root.exe ${cppfilename} -b -q )
   set_tests_properties(${cppfname}CPP PROPERTIES TIMEOUT 5000 ) 
endforeach() 

file(GLOB cxxfiles "${CMAKE_INSTALL_PREFIX}/examples/*/*.cxx")
foreach(cxxfilename ${cxxfiles})
   get_filename_component(cxxfname ${cxxfilename} NAME_WE)
   string(CONCAT cxxfname ostap.examples_cxx_ ${cxxfname})
   add_test (NAME ${cxxfname}CXX WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary2 COMMAND root.exe ${cxxfilename} -b -q )
   set_tests_properties(${cxxfname}CXX PROPERTIES TIMEOUT 5000 ) 
endforeach()


