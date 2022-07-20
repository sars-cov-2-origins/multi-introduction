A specific version of BEAST is required to run these XMLs. Please follow the below instructions to install and use the appropriate version. 

1. Install the latest versions of BEAST and BEAGLE, following the instructions from http://beast.community/installing. 

2. Install Apache Ant, following the instructions from https://ant.apache.org/manual/install.html.

3. Download the BEAST `sars-cov-2-origins` release from https://github.com/beast-dev/beast-mcmc/releases/sars-cov-2-origins.

4. Move to the top of the directory for the BEAST `sars-cov-2-origins` release.

5. Execute: `ant` (this will build a JAR file called `./build/dist/beast.jar`).

6. Copy `beast.jar` from the sars-cov-2-origins release on top of the `beast.jar` file in the `lib` subdirectory of the BEAST installation from step 1. 

7. Run BEAST from step 1 as you normally would.
