set(APP_NAME "dynamics.${STORMM_APP_SUFFIX}")

add_executable(${APP_NAME}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/setup.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/setup.h)

target_sources(${APP_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/simulator.cpp)
      
target_link_libraries(${APP_NAME} ${PROJECT_NAME})

install(TARGETS ${APP_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
