This is the code and data for the model-subject comparisons in the paper [There's Waldo! A Normalization Model of Visual Search Predicts Single-Trial Human Fixations in an Object Search Task](https://academic.oup.com/cercor/article/26/7/3064/1745172/There-s-Waldo-A-Normalization-Model-of-Visual) (Miconi, Groomes and Kreiman, Cereb Cortex (2015) 26 (7): 3064-3082). 

Directory `subjectsdata` contains the human subject data from the visual search task. Be sure to read the instructions in this directory. In particular, note that several files contain incomplete or buggy data!

To run the model:

- Compile `hmaxC2b.cpp` (e.g. with `g++ -O2 -Wall hmaxC2b.cpp -o hmaxC2b`)

- In MATLAB, run `makepics'. This will generate the pictures on which the model will be applied (these images are identical to the ones seen by human subjects, except for overall size).

- Run `./hmaxC2b`. this will appl ythe HMAX code to the images, which will generate S2b and C2b outputs in the `out` directory. Note that this will take a lot of time (on the order of a day) and generate a lot of data (~60GB).

- In MATLAB, run `attentionMapC`. This will apply the visual search model to the HMAX outputs. The output will be stored in `modeloutput.mat`.

- In MATLAB, run `analysetest.m`, which will generate (some of) the figures from the paper.

