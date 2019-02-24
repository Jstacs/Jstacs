# Description

Sequence analysis is one of the major subjects of bioinformatics. Several existing libraries combine the representation of biological sequences with exact and approximate pattern matching as well as alignment algorithms. We present Jstacs, an open source Java library, which focuses on the statistical analysis of biological sequences instead. Jstacs comprises an efficient representation of sequence data and provides implementations of many statistical models with generative and discriminative approaches for parameter learning. Using Jstacs, classifiers can be assessed and compared on test datasets or by cross-validation experiments evaluating several performance measures. Due to its strictly object-oriented design Jstacs is easy to use and readily extensible.

For more information including an API documentation, code examples, FAQs, binaries, and a cookbook visit http://www.jstacs.de.

# Organization of the library

Jstacs core classes may be found in sub-packages of de.jstacs.
Individual projects using and extending these core classes for specific applications are located in packages projects.

A list of projects that are based on Jstacs, including binaries documentation of user parameters is available at http://jstacs.de/index.php/Projects.

Building upon Jstacs, [JstacsFX](https://github.com/Jstacs/JstacsFX) visualizes parameters and results in a JavaFX-based GUI that is built upon the generic de.jstacs.tools.JstacsTool class.

# Licensing information

Jstacs is free software: you can redistribute it and/or modify under the terms of the GNU General Public License version 3 or (at your option) any later version as published by the Free Software Foundation.

For more information, please read COPYING.txt.

