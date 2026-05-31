using TestItemRunner

@run_package_tests filter=ti -> !(manualsuite in ti.tags)