Description
-----------
Sequence analysis is one of the major subjects of bioinformatics. Several existing libraries combine the representation of biological sequences with exact and approximate pattern matching as well as alignment algorithms. We present Jstacs, an open source Java library, which focuses on the statistical analysis of biological sequences instead. Jstacs comprises an efficient representation of sequence data and provides implementations of many statistical models with generative and discriminative approaches for parameter learning. Using Jstacs, classifiers can be assessed and compared on test datasets or by cross-validation experiments evaluating several performance measures. Due to its strictly object-oriented design Jstacs is easy to use and readily extensible.

For more information including an API, code examples, FAQs, and a forum visit http://www.jstacs.de.


Licensing Information
--------------------- 
Jstacs is free software: you can redistribute it and/or modify under the terms of the GNU General Public License version 3 or (at your option) any later version as published by the Free Software Foundation.

For more information please read COPYING.txt. 


Libraries
---------
The development of Jstacs greatly profited from the following libraries 
* BioJava, 
* numerical methods from JTEM,
* Rserve, and
* RandomNumberGenerator was adapted from Sundar Dorai-Raj. 


Getting started
---------------
If you like to use the Jstacs binaries in your own project, you must include jstacs-<version>.jar and all JARs in the lib-directory into your classpath, e.g.

java -cp ./:./jstacs-<version>.jar:./lib/biojava-live.jar:./lib/bytecode.jar:./lib/numericalMethods.jar:./RClient.jar <test-class>

If you want to use the class UserTime, which uses native code, you must also set the library-path to the directory where the dynamic library resides, e.g.

java -Djava.library.path=./native/ -cp <see above> <test-class>


Bug reports & Feature requests
------------------------------ 
You can submit bug reports and feature requests via the Jstacs trac at http://www.jstacs.de. Before you open a new bug ticket, please check if that bug has already been submitted in the list of existing tickets.
In the Jstacs trac, we also provide a forum for discussions about Jstacs. 