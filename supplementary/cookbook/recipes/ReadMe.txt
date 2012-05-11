
= Introduction =

These code examples shall show how you can use Jstacs to get several tasks done.
Each of these code examples must be compiled with Jstacs and the required libraries in the classpath, e.g.,

 javac -cp jstacs-2.0.jar AlphabetCreation.java

or

 javac -cp jstacs-2.0.jar:RClient-0.6.7.jar RserveTest.java

If the libraries are not located in the same directory as the Java source file, you need to adapt the paths to the libraries Jars.

The code examples
* AlphabetCreation,
* TrainPWM,
* TrainHomogeneousMM,
* CreateMixtureModel,
* AnalyseDataWithDifferentModels,
* DeNovoSunflower,
* TrainClassifier,
* CreateMSPClassifier,
* CrossValidation,
* CurvePlotter,
* DataLoader,
* GenDisMixClassifierTest,
* TrainSMBasedClassifierTest, and
* RserveTest
contain a main-method. Hence, after compilation, you can directly run them as a program, e.g.,

 java -cp .:jstacs-2.0.jar:RClient-0.6.7.jar RserveTest

The remaining two examples
* HomogeneousMarkovModel and
* PositionWeightMatrixDiffSM
show how to implement new TrainableStatisticalModels and DifferentiableStatisticalModels and do not contain a
main-method.

= Examples =

== AlphabetCreation ==
In this example, we create a new ComplementableDiscreteAlphabet using the generic implementation. 
We then use this Alphabet to create a Sequence and compute its reverse complement.

=== Compile ===
 javac -cp jstacs-2.0.jar AlphabetCreation.java

=== Run ===
 java -cp .:jstacs-2.0.jar AlphabetCreation

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== TrainPWM ==
In this example, we show how to load sequence data into Jstacs and how to learn a position weight 
matrix (inhomogeneous Markov model of order 0) on these data.

=== Compile ===
javac -cp jstacs-2.0.jar TrainPWM.java

=== Run ===
java -cp .:jstacs-2.0.jar TrainPWM fg.fa
where "fg.fa" may be replaced by any FastA file containing fixed-length sequences.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== TrainHomogeneousMM ==
In this example, we show how to load sequence data into Jstacs and how to learn a 
homogeneous Markov model of order 1 on these data.

=== Compile === 
javac -cp jstacs-2.0.jar TrainHomogeneousMM.java

=== Run ===
java -cp .:jstacs-2.0.jar TrainHomogeneousMM fg.fa
where "fg.fa" may be replaced by any FastA file.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== GenerateData ==
In this example, we show how to learn a homogeneous Markov model of order 2 
from data (similar to the previous example), and use the learned model to generate 
new data following the same distribution as the original data.

=== Compile === 
javac -cp jstacs-2.0.jar GenerateData.java

=== Run ===
java -cp .:jstacs-2.0.jar GenerateData fg.fa generated.txt
where "fg.fa" may be replaced by any FastA file, and after running the program,
"generated.txt" contains the generated data as plain text.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== CreateMixtureModel ==
In this example, we show how to load sequence data into Jstacs and how to learn a 
mixture model of two position weight matrices on these data using the expectation 
maximization algorithm.

=== Compile === 
javac -cp jstacs-2.0.jar CreateMixtureModel.java

=== Run ===
java -cp .:jstacs-2.0.jar CreateMixtureModel fg.fa
where "fg.fa" may be replaced by any FastA file containing fixed-length sequences.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== AnalyseDataWithDifferentModels ==
In this example, we show how to use the TrainableStatisticalModelFactory to create 
inhomogeneous and homogeneous Markov models, and Bayesian trees, and how to learn 
these models on a common data set.

=== Compile === 
javac -cp jstacs-2.0.jar AnalyseDataWithDifferentModels.java

=== Run ===
java -cp .:jstacs-2.0.jar AnalyseDataWithDifferentModels fg.fa
where "fg.fa" may be replaced by any FastA file containing fixed-length sequences.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== DeNovoSunflower ==
In this example, we show how to use the HMMFactory to create a sunflower hidden 
Markov model (HMM) with two motifs of different lengths. We show how to train the 
sunflower HMM on input data, which are typically long sequences containing an 
over-represented motif.
After training the HMM, we show how to compute and output the Viterbi paths for all 
sequences, which give an indication of the position of motif occurrences.

=== Compile === 
javac -cp jstacs-2.0.jar DeNovoSunflower.java

=== Run ===
java -cp .:jstacs-2.0.jar DeNovoSunflower promoters.fa
where "promoters.fa" may be replaced by any FastA file containing sequences with
a hidden motif. The example data set "promoters.fa" contains artificial sequences
drawn from a uniform distribution with hidden binding sites extracted from Jaspar
(http://jaspar.genereg.net/).

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== TrainClassifier ==
In this example, we show how to train a classifier based on a position weight matrix 
model and a homogeneous Markov model on training data, and how to use the trained 
classifier to classify sequences.

=== Compile === 
javac -cp jstacs-2.0.jar TrainClassifier.java

=== Run ===
java -cp .:jstacs-2.0.jar TrainClassifier fg.fa bg.fa
where "fg.fa" may be replaced by any FastA file containing fixed-length sequences,
and "bg.fa" may be replaced by any FastA file.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== CreateMSPClassifier ==

In this example, we show how to use the \DiffSMFactory~to create a position weight 
matrix and how to learn a classifier based on two position weight matrices using the 
discriminative maximum supervised posterior principle.

=== Compile === 
javac -cp jstacs-2.0.jar CreateMSPClassifier.java

=== Run ===
java -cp .:jstacs-2.0.jar CreateMSPClassifier fg.fa bg.fa
where "fg.fa" may be replaced by any FastA file containing fixed-length sequences,
and "bg.fa" may be replaced by any FastA file.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== DataLoader ==
In this example, we show different ways of creating a DataSet in Jstacs from plain text and FastA 
files and using the adaptor to BioJava.

=== Compile ===
 javac -cp jstacs-2.0.jar:biojava-live.jar:bytecode.jar DataLoader.java

=== Run ===
 java -cp .:jstacs-2.0.jar:biojava-live.jar:bytecode.jar DataLoader ./
where the files "myfile.fa" and "myfile.txt" must be located in the working directory.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== TrainSMBasedClassifierTest ==
In this example, we show how to create a TrainSMBasedClassifier using to position weight matrices, 
train this classifier, classify previously unlabeled data, store the classifier to its XML representation, 
and load it back into Jstacs.

=== Compile ===
 javac -cp jstacs-2.0.jar TrainSMBasedClassifierTest.java

=== Run ===
 java -cp .:jstacs-2.0.jar TrainSMBasedClassifierTest ./ fg.fa bg.fa unknown.fa
where the files "fg.fa", "bg.fa", and "unknown.fa" must be located in the working directory.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== GenDisMixClassifierTest ==
In this example, we show how to create GenDisMixClassifiers using two position weight matrices. 
We show how GenDisMixClassifiers can be created for all basic learning principles (ML, MAP, MCL, MSP), 
and how these classifiers can be trained and assessed.

=== Compile ===
 javac -cp jstacs-2.0.jar:numericalMethods.jar GenDisMixClassifierTest.java

=== Run ===
 java -cp .:jstacs-2.0.jar:numericalMethods.jar GenDisMixClassifierTest fg.fa bg.fa

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== RserveTest ==
Here, we show a number of examples how R can be used from within Jstacs using RServe.

=== Compile ===
 javac -cp jstacs-2.0.jar:RClient-0.6.7.jar RserveTest.java

=== Run ===
 java -cp .:jstacs-2.0.jar:RClient-0.6.7.jar RserveTest ./
where plots are created in the working directory. For this program, 
a running Rserve on the local computer with standard port is required.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== CurvePlotter ==
In this example, we show how a classifier (loaded from disk) can be assessed on test data, 
and how we can plot ROC and PR curves of this classifier and test data set.

=== Compile ===
 javac -cp jstacs-2.0.jar:RClient-0.6.7.jar CurvePlotter.java

=== Run ===
 java -cp .:jstacs-2.0.jar:RClient-0.6.7.jar CurvePlotter ./ fg.fa bg.fa
where the files "fg.fa", "bg.fa", and "myClassifier.xml" must be located in the working directory. 
The classifier in "myClassifier.xml" can be created by running the code example "TrainSMBasedClassifierTest".
For this program, a running Rserve on the local computer with standard port is required.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== CrossValidation ==
In this example, we show how we can compare classifiers built on different types of models and using different 
learning principles in a cross validation. Specifically, we create a position weight matrix, use that matrix to 
create a mixture model, and we create an inhomogeneous Markov model of order 3. 
We do so in the world of TrainableStatisticalModels and in the world of DifferentiableStatisticalModels. 
We then use the mixture model as foreground model and the inhomogeneous Markov model as the background model
when building classifiers. The classifiers are learned by the generative MAP principle and the discriminative 
MSP principle, respectively. 
We then assess these classifiers in a 10-fold cross validation.

=== Compile ===
 javac -cp jstacs-2.0.jar:numericalMethods.jar CrossValidation.java

=== Run ===
 java -cp .:jstacs-2.0.jar:numericalMethods.jar CrossValidation fg.fa bg.fa
where "fg.fa" may be replaced by any FastA file containing fixed-length sequences,
and "bg.fa" may be replaced by any FastA file.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

== HomogeneousMarkovModel ==
In this example, we show how to implement a new TrainableStatisticalModel. Here, we implement a simple homogeneous 
Markov models of order 0 to focus on the technical side of the implementation.

=== Compile ===
 javac -cp jstacs-2.0.jar HomogeneousMarkovModel.java
 
== PositionWeightMatrixDiffSM ==
In this example, we show how to implement a new DifferentiableStatisticalModel. Here, we implement a simple position 
weight matrix, i.e., an inhomogeneous Markov model of order 0. Since we want to use this position weight matrix in 
numerical optimization, we parameterize it in the so called "natural parameterization".
Since we use a product-Dirichlet prior on the parameters, we transformed this prior to the parameterization we use.

=== Compile ===
 javac -cp jstacs-2.0.jar PositionWeightMatrixDiffSM.java