set(kritahalftone_SOURCES
    KisHalftoneFilter.cpp
    KisHalftoneFilterConfiguration.cpp
    KisHalftoneConfigWidget.cpp
    KisHalftoneConfigPageWidget.cpp
)

ki18n_wrap_ui(kritahalftone_SOURCES
    KisHalftoneConfigWidget.ui
    KisHalftoneConfigPageWidget.ui
)

kis_add_library(kritahalftone MODULE ${kritahalftone_SOURCES})
target_link_libraries(kritahalftone kritaui)
install(TARGETS kritahalftone DESTINATION ${KRITA_PLUGIN_INSTALL_DIR})
install( FILES
    halftone.action
DESTINATION  ${KDE_INSTALL_DATADIR}/krita/actions)
