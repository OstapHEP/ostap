
install ( FILES       ostap 
          DESTINATION scripts
          COMPONENT   scripts 
          PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install ( FILES       check-dependencies 
          DESTINATION scripts  
          COMPONENT   scripts 
          PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install ( FILES       sync-dirs 
          DESTINATION scripts  
          COMPONENT   scripts 
          PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install ( FILES       pplaunch
          DESTINATION scripts
          COMPONENT   scripts 
          PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/thisostap.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/thisostap.sh" @ONLY
  )
    
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/.ostaprc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/.ostaprc"     @ONLY
  )

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/rootlogon.C.in"
  "${CMAKE_CURRENT_BINARY_DIR}/.rootlogon.C" @ONLY
  )

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/rootlogon.C.in"
  "${CMAKE_CURRENT_BINARY_DIR}/../examples/.rootlogon.C" @ONLY
  )
    
file    ( TOUCH "${CMAKE_CURRENT_SOURCE_DIR}/.footprints" ) 
        
install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/.rootlogon.C  
          DESTINATION .
          COMPONENT   scripts 
          PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install( FILES ${CMAKE_CURRENT_BINARY_DIR}/thisostap.sh  
         DESTINATION . 
         COMPONENT   scripts 
         PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install( FILES ${CMAKE_CURRENT_BINARY_DIR}/.ostaprc   
         DESTINATION . 
         COMPONENT   scripts 
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ )

install ( FILES ${CMAKE_CURRENT_SOURCE_DIR}/.footprints   
          DESTINATION .  
          COMPONENT   scripts 
          PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_WRITE WORLD_READ )

##  senseless? 
set( $ENV{PATH} ${CMAKE_INSTALL_PREFIX}/bin:$ENV{PATH} ) 
