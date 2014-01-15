install(CODE "MESSAGE(
\"----------------------------------------------------------------------------
TerraFERMA has now been installed in

    ${CMAKE_INSTALL_PREFIX}

Before using please update your environment variables. 

This can be done easily using the helper file 'terraferma.conf' which sets 
the appropriate variables (for users of the Bash shell).

To update your environment variables, run the following command:

    source ${CMAKE_INSTALL_PREFIX}/share/terraferma/terraferma.conf

For future reference, we recommend that you add this command to your
configuration (.bashrc, .profile or similar).

Alternatively, an environment module has been installed in:

    ${CMAKE_INSTALL_PREFIX}/share/terraferma/terraferma.configmodule

Either copy it to an appropriate location in, or add that location to, your 
\${MODULESPATH} before running 'module load terraferma.configmodule'.
----------------------------------------------------------------------------\")")

