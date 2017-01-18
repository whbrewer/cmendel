MENDEL'S ACCOUNTANT 

Mendel's Accountant (MENDEL) is an advanced numerical simulation program for 
modeling genetic change over time and was developed collaboratively by Sanford,
Baumgardner, Brewer, Gibson and ReMine.

This C-version of the Mendel engine was developed in 2008 by Wesley Brewer 
and Huang Jie with the grateful support of the FMS Foundation, especially
Dr. John Sanford at Cornell University.

To compile the serial version:

> make serial

In order to make the Parallel Version, mpich2 must first be installed, then:

> make 

To run Mendel's Accountant, modify the file mendel.in to set the parameters
in which you are interested.  Then, type:

> ./mendel

--
"...the truth shall make you free..."

