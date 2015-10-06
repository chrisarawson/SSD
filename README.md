# SSD
Code to return an SSD with protection levels.

Code returns:
<ul>
<li> <fit> An object with the parameters of the lognormal fit to the data</li>
<li> <hcs> a vector of the 90th, 95th and 99th percent values from the lognormal fit</li>
<li> <protValTab> A table of the mean and 95% CIS for the 99th, 95th and 90th percent values based on bootstrapped values. These are likely to be different from the values from <fit></li>
<li> A plot of the data points so the user can check for doubles etc</li>
<li> A plot of the data broken by a factor, the bootstrapped lognormal fit to the data and 95% confidence intervals.</li>
 

Sample data is contained in the file LC50.csv

Written by C.Rawson based on code from Eduard Szocs (http://edild.github.io/ssd/)
