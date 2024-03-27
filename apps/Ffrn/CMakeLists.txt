set(APP_NAME "ffrefine.${STORMM_APP_SUFFIX}")

add_executable(${APP_NAME}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ffrefine.cpp)

target_link_libraries(${APP_NAME} ${PROJECT_NAME})

install(TARGETS ${APP_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
