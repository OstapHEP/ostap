
install(FILES ostap 
        DESTINATION scripts
        PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(FILES ostap-check-dependencies 
        DESTINATION scripts  
        PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(FILES pplaunch
        DESTINATION scripts
        PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/thisostap.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/thisostap.sh" @ONLY
  )
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/.ostaprc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/.ostaprc" @ONLY
  )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/thisostap.sh  
        DESTINATION . 
        PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/.ostaprc   
        DESTINATION . 
        PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

##  senseless? 
set( $ENV{PATH} ${CMAKE_INSTALL_PREFIX}/bin:$ENV{PATH} ) 

