package supplementary.codeExamples;

import java.io.PrintWriter;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifier.AbstractScoreBasedClassifier;
import de.jstacs.classifier.ConfusionMatrix;
import de.jstacs.classifier.MeasureParameters;
import de.jstacs.classifier.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifier.MeasureParameters.Measure;
import de.jstacs.classifier.assessment.RepeatedHoldOutAssessParameterSet;
import de.jstacs.classifier.assessment.RepeatedHoldOutExperiment;
import de.jstacs.classifier.modelBased.ModelBasedClassifier;
import de.jstacs.classifier.scoringFunctionBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifier.scoringFunctionBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifier.scoringFunctionBased.logPrior.CompositeLogPrior;
import de.jstacs.classifier.scoringFunctionBased.msp.MSPClassifier;
import de.jstacs.data.DNASample;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.models.Model;
import de.jstacs.models.VariableLengthWrapperModel;
import de.jstacs.models.discrete.inhomogeneous.BayesianNetworkModel;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.models.discrete.inhomogeneous.parameters.BayesianNetworkModelParameterSet;
import de.jstacs.results.ListResult;
import de.jstacs.results.ResultSet;
import de.jstacs.scoringFunctions.directedGraphicalModels.BayesianNetworkScoringFunction;
import de.jstacs.scoringFunctions.directedGraphicalModels.BayesianNetworkScoringFunctionParameterSet;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.utils.REnvironment;


/**
 * This class implements a main that shows some features of Jstacs including models for generative and discriminative learning,
 * creation of classifiers, evaluation of classifiers, hold-out sampling, and binding site prediction.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class MiMBExample {

	/**
	 * @param args only the first parameter will be used; it determines the home directory
	 */
	public static void main( String[] args ) throws Exception {
	
		String home = args[0]+System.getProperty( "file.separator" );
		
		/* read data */
		
		//read foreground and background data set form FastA-files
		Sample fgData = new DNASample(home+"foreground.fa");
		Sample bgData = new DNASample(home+"background.fa");
		
		/* generative part */
		
		//create set of parameters for foreground model
		BayesianNetworkModelParameterSet pars = new BayesianNetworkModelParameterSet(
				fgData.getAlphabetContainer(),//used alphabets
				fgData.getElementLength(),//element length == sequence length of each sequence in the sample
				4,//ESS == equivalent sample size (has to be non-negative)
				"fg model",//user description of the model
				ModelType.IMM,//type of statistical model, here an inhomogeneous Markov model (IMM)
				(byte)0,// model order, here 0, so we get an IMM(0) == PWM = position weight matrix
				LearningType.ML_OR_MAP//how to learn the parameters, depends on ESS; for ESS=0 it is ML otherwise MAP
		);
		//create foreground model from these parameters
		Model fgModel = new BayesianNetworkModel(pars);
		
		//analogously, create the background model
		BayesianNetworkModelParameterSet pars2 = new BayesianNetworkModelParameterSet(
				fgData.getAlphabetContainer(),
				fgData.getElementLength(),
				1024,
				"bg model",
				ModelType.IMM,
				(byte)0,
				LearningType.ML_OR_MAP
		);		
		Model bgModel = new BayesianNetworkModel(pars2);
		bgModel = new VariableLengthWrapperModel(bgModel);
		
		//create generative classifier from the models defined before
		ModelBasedClassifier cl = new ModelBasedClassifier(fgModel, bgModel);	
		
		cl.train( fgData, bgData );
		
		/* discriminative part */
		
		//create set of parameters for foreground scoring function
		BayesianNetworkScoringFunctionParameterSet parsD = new BayesianNetworkScoringFunctionParameterSet(
				fgData.getAlphabetContainer(),//used alphabets
				fgData.getElementLength(),//element length == sequence length of each sequence in the sample
				4,//ESS == equivalent sample size (has to be non-negative)
				true,//use MAP-parameters as start values of the optimization
				new InhomogeneousMarkov(0) //the statistical model, here an IMM(0) == PWM
		);
		//create foreground scoring function from these parameters
		BayesianNetworkScoringFunction fgFun = new BayesianNetworkScoringFunction(parsD);
		
		//analogously, create the background scoring function
		BayesianNetworkScoringFunctionParameterSet parsDbg = new BayesianNetworkScoringFunctionParameterSet(
				fgData.getAlphabetContainer(),
				fgData.getElementLength(),
				1024,
				true,
				new InhomogeneousMarkov(0));
		BayesianNetworkScoringFunction bgFun = new BayesianNetworkScoringFunction(parsDbg);
		
		//create set of parameter for the discriminative classifier
		GenDisMixClassifierParameterSet clPars = new GenDisMixClassifierParameterSet(
				fgData.getAlphabetContainer(),//used alphabets
				fgData.getElementLength(),//element length == sequence length of each sequence in the samples
				Optimizer.QUASI_NEWTON_BFGS,//determines the algorithm for numerical optimization
				1E-6,//epsilon to stop the numerical optimization
				1E-6,//epsilon to stop the line search within the numerical optimization
				1,//start step width in the numerical optimization 
				false,//a switch that decides whether to use only the free parameters or all parameters
				KindOfParameter.PLUGIN,//a switch to decide which start parameters to choose
				true,//a switch that states the objective function will be normalized
				1//number of threads used during optimization
		);
		//create discriminative classifier from the parameters and the scoring function defined before
		MSPClassifier cll = new MSPClassifier( 
				clPars,//the parameters of the classifier
				new CompositeLogPrior(),//the used prior, to obtain MCL use null
				fgFun,//the scoring function for the foreground class
				bgFun//the scoring function for the background class
		);
		
		//train the discriminative classifier
		cll.train( fgData, bgData );
		
		/* performance measures */
		
		//partition data
		Sample[] fgSplit = bisect( fgData, fgData.getElementLength() );
		Sample[] bgSplit = bisect( bgData, fgData.getElementLength() );

		Sample fgTest = fgSplit[1];
		Sample bgTest = bgSplit[1];
		
		//train the generative classifier
		cl.train( fgSplit[0], bgSplit[0] );
		
		//fill a confusion matrix
		ConfusionMatrix confMatrix = cl.test( fgTest, bgTest );
		
		//read the entries of the table
		double tp = confMatrix.getCountsFor(0, 0);
		double fn = confMatrix.getCountsFor(1, 0);
		double tn = confMatrix.getCountsFor(1, 1);
		double fp = confMatrix.getCountsFor(0, 1);
		
		double p = tp+fn;
		double barp = tp+fp;
		double n = tn+fp;
		double barn = tn+fn;
		
		System.out.println("TP = "+tp+"\t\tFP = "+fp+"\t\tbarp = "+barp+"\n" +
				"FN = "+fn+"\t\tTN = "+tn+"\t\tbarn = "+barn+"\n" +
				"p = "+p+"\t\tn = "+n+"\t\tN' = "+(n+p));
		
		//compute the measures
		double sn = tp/p;
		double ppv = tp/barp;
		double fpr = fp/n;
		double sp = tn/n;
		double cr = (tp+tn)/(n+p);
		
		System.out.println("cr = "+cr+"\nSn = "+sn+"\nppv = "+ppv+"\nSp = "+sp+"\nfpr = "+fpr+"\n");
		
		//define the measures that shall be evaluated
		MeasureParameters mp = new MeasureParameters(
				true,//evaluate all performance measures
				0.999,//use specificity of 0.999 to measure the sensitivity
				0.95,//use sensitivity of 0.95 to measure the false positive rate 
				0.95//use sensitivity of 0.95 to measure the positive predictive value
		);
		
		//evaluates the classifier
		ResultSet rs = cl.evaluateAll(
				mp,//defines the measures that will be evaluated
				true,//allows to throw an exception if a measure can not be computed
				fgTest,//the test data for the foreground class
				bgTest//the test data for the background class
		);
		System.out.println(rs); 

		//plot ROC and PR curve
		DoubleTableResult roc = (DoubleTableResult)rs.getResultAt( rs.findColumn( Measure.ReceiverOperatingCharacteristicCurve.getNameString() ) );
		DoubleTableResult pr = (DoubleTableResult)rs.getResultAt( rs.findColumn( Measure.PrecisionRecallCurve.getNameString() ) );
		
		REnvironment re = null;
		//you need to have a Rserve running
		try {
			re = new REnvironment(
					"localhost",//server name
					"",//user name
					""//password
				);
			
			String snfpr = "points( " + fpr + ", " + sn + ", col=" + 1 	+ ", pch=" +4 + ", cex=2, lwd=3 );\n";
			String ppvsn = "points( " + sn + ", " + ppv + ", col=" + 1 	+ ", pch=" +4 + ", cex=2, lwd=3 );\n";
			re.voidEval( "p<-palette();p[8]<-\"gray66\";palette(p);" );
			
			re.plotToPDF( DoubleTableResult.getPlotCommands( re, null, new int[]{8}, roc ).toString()+"\n"+snfpr,4,4.5, home+"roc.pdf",true);
			re.plotToPDF( DoubleTableResult.getPlotCommands( re, null, new int[]{8}, pr ).toString()+"\n"+ppvsn, 4,4.5, home+"pr.pdf",true);
		} catch( Exception e ) {
			System.out.println( "could not plot the curves" );
		} finally {
			if( re != null ) {
				re.close();
			}
		}
		separator();
		
		/* hold-out sampling */
		
		//define the measures that shall be evaluated
		mp = new MeasureParameters(
				false,//only evaluate numerical performance measures
				0.999,//use specificity of 0.999 to measure the sensitivity
				0.95,//use sensitivity of 0.95 to measure the false positive rate 
				0.95//use sensitivity of 0.95 to measure the positive predictive value
		);
		
		//create the parameters for the hold-out sampling
		RepeatedHoldOutAssessParameterSet parsA = new RepeatedHoldOutAssessParameterSet(
			Sample.PartitionMethod.PARTITION_BY_NUMBER_OF_SYMBOLS,//defines the way of splitting the data
			fgData.getElementLength(), //defines the length of the elements in the test data set
			true,//a switch that decides whether to throw an exception if a performance measure can not be evaluated
			1000,//the number of samplings
			new double[]{0.1,0.1}//the partition of the data of each class the will be used for testing
		);
		
		//creates an hold-out experiment for two classifiers (the generative and the discriminative classifier)
		RepeatedHoldOutExperiment exp = new RepeatedHoldOutExperiment(cl,cll); 
		//does the experiment a stores the results in a ListResult
		ListResult lr = exp.assess(
				mp,//the measures that will be computed
				parsA,//the parameters for the experiment
				fgData,//the foreground data
				bgData//the background data
		);
		
		System.out.println(lr);
		separator();
		
		/* prediction */
		
		//re-train discriminative classifier
		cll.train( fgData, bgData );
		
		//load data for prediction
		Sample promoters = new DNASample(home+"human_promoters.fa");
		
		//find best possible binding site
		int si=0,id=0;
		double llr, max=Double.NEGATIVE_INFINITY;
		
		PrintWriter out = new PrintWriter( home+"/allscores.txt" );
		
		int i = 0;
		//check all sequence
		for( Sequence seq : promoters ){
			//check each possible start position
			for(int l=0;l<seq.getLength()-cll.getLength()+1;l++){
				Sequence sub = seq.getSubSequence( l, cll.getLength() );
				//compute likelihood ratio
				llr = cll.getScore( sub, 0 ) - cll.getScore( sub, 1 );
				out.print( llr + "\t" );
				if(llr > max){
					//set new best likelihood ratio, sequence, and site
					max = llr;
					si = i;
					id = l;
				}
			}
			out.println();
			i++;
		}
		out.close();
		
		//write a file for the best prediction
		Sequence bestSequence = promoters.getElementAt( si );
		out = new PrintWriter(home+"/scores.txt");
		//write the sequence
		out.println(bestSequence.toString("\t", id-30,id+30));
		
		//write the log likelihood ratio that can be used to plot a profile
		for(int l=id-30;l<id+30;l++){
			Sequence site = bestSequence.getSubSequence( l, cll.getLength() );
			out.print((cll.getScore( site, 0 ) - cll.getScore( site, 1 ))+"\t");
		}
		out.println();
		out.close();
	}
	
	//method for obtaining reproduceably the same split of some given data
	private static Sample[] bisect(Sample data, int l) throws Exception {
		int mid = data.getNumberOfElements()/2;
		return new Sample[]{
				getSubSample(data, 0, mid, "train",l),
				getSubSample(data, mid, data.getNumberOfElements(), "test",l)
		};
	}
	
	//creates a sample from a specific part of the data
	private static Sample getSubSample( Sample data, int start, int end, String annotation, int l ) throws Exception {
		//copy the sequences into an array
		Sequence[] seqs = new Sequence[end-start];
		for(int i=0;i<seqs.length;i++){
			seqs[i] = data.getElementAt( i+start );
		}
		return new Sample( new Sample( annotation, seqs ), l );
	}
	
	//prints a separator
	private static void separator() {
		for( int i = 0; i < 50; i++) {
			System.out.print("=");
		}
		System.out.println();
	}
}