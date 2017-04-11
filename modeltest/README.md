This is the code for the tests of model performance.


To use this code, you will need MATLAB and a C++ compiler.

How to compile and run the code:

- With MATLAB, run `makepics.m`. This will generate the pictures on which model performance is tested. The pictures will be stored in the `images` directory.

- Compile `hmaxC2b.cpp`, which implements the HMAX model (e.g. `g++ -Wall hmaxC2b.cpp -o hmaxC2b`).

- Run `./hmaxC2b`. This will apply the HMAX model to all the (raw) images listed in the file `listinputfiles.txt` (which should be identical to the results of `find images/*.raw` *after* you run `makepics.m`). This will produce S2b and C2b value files (see paper). This will take a LOT of time (perhaps a day) and generate a LOT of data (on the order of 59 GB).

- With MATLAB, run `runModelArgs('array')` (which will apply the model to Array-type images) and then 'runModelArgs('unpred')` (which will apply the model to natural-background images). The result of the model is the text output produced by this program, which will simply report success rate for each successive fixation.


