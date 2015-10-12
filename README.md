# SSD
Code to return an SSD with protection levels.

Code returns:
<ul>
<li> <b>fit</b> An object with the parameters of the lognormal fit to the data</li>
<li> <b>hcs</b> a vector of the 90th, 95th and 99th percent values from the lognormal fit</li>
<li> <b>protValTab</b> A table of the mean and 95% CIS for the 99th, 95th and 90th percent values based on bootstrapped values. These are likely to be different from the values from <b>fit</b></li>
<li> A plot of the data points so the user can check for doubles etc</li>
<li> A plot of the data broken by a factor, the bootstrapped lognormal fit to the data and 95% confidence intervals.</li>
 
Also generated is an object <b>mySumm</b> that summarises the output.
Calling >mySumm will output the summary.

Sample data is contained in the file LC50.csv collated by Amelia Wenger

Written by C.Rawson based on code from Eduard Szocs (http://edild.github.io/ssd/)
