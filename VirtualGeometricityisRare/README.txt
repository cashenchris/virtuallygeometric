The scripts in this folder can be used to recreate the graphs from the paper "Virtual Geometricity is Rare" by Christopher H. Cashen and Jason F. Manning.

Syntax is:
./*experiment.py rank minwordlength maxwordlength numberoftrials

For example, in rank 2, doing:

./VGexperiment.py 2 1 15 100
./fullwordexperiment.py 2 1 100 500
./extract_and_plot2.py 2 rank2gvg.txt rank2fullwords.txt

runs experiments and creates graphs for rank 2 and saves them as rank2gvgfull.pdf and rank2loggvg.pdf. 

Exponential approximation only takes into account the data when the proportion of geometric/non-full words is <=.5. Lack of such data results in:
"RuntimeError: data does not go below .5"
Try another experiment with larger maxwordlength.

Exponential approximation also fails with error:
"ValueError: array must not contain infs or NaNs" 
if proportion of geometric/non-full words for some length is zero. 
Try another experiment with larger numberoftrials. 
