package supplementary.codeExamples;

import java.io.OutputStream;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Random;

import org.biojava.bio.seq.SequenceIterator;
import org.biojavax.bio.db.RichSequenceDB;
import org.biojavax.bio.db.ncbi.GenbankRichSequenceDB;
import org.biojavax.bio.seq.RichSequenceIterator;

import de.jstacs.DataType;
import de.jstacs.Storable;
import de.jstacs.algorithms.alignment.Alignment;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.cost.Costs;
import de.jstacs.algorithms.alignment.cost.SimpleCosts;
import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.DifferentiableFunction;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.NumericalDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.CombinedCondition;
import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.classifiers.assessment.ClassifierAssessment;
import de.jstacs.classifiers.assessment.KFoldCrossValidation;
import de.jstacs.classifiers.assessment.KFoldCrossValidationAssessParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure;
import de.jstacs.classifiers.performanceMeasures.AucPR;
import de.jstacs.classifiers.performanceMeasures.AucROC;
import de.jstacs.classifiers.performanceMeasures.NumericalPerformanceMeasureParameterSet;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet;
import de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.data.alphabets.Alphabet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.GenericComplementableDiscreteAlphabet;
import de.jstacs.data.bioJava.BioJavaAdapter;
import de.jstacs.data.sequences.PermutedSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.annotation.MotifAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SequenceScoringParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.MeanResultSet;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.sampling.VarianceRatioBurnInTest;
import de.jstacs.sampling.VarianceRatioBurnInTestParameterSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.IndependentProductDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.btMeasures.BTExplainingAwayResidual;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.MixtureDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM.InitMethod;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.motif.ExtendedZOOPSDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.DifferentiableStatisticalModelWrapperTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.homogeneous.HomogeneousMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.homogeneous.parameters.HomMMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.BayesianNetworkTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.BayesianNetworkTrainSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.HMMFactory;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.DifferentiableHigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.HigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete.DiscreteEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.BaumWelchParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.HMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.NumericalHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.MixtureTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.StrandTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.motif.ZOOPSTrainSM;
import de.jstacs.utils.REnvironment;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.SubclassFinder;
import de.jstacs.utils.ToolBox;


public class NewCodeExampleTest {

	/**
	 * @param args
	 */
	public static void main( String[] args ) {
		// TODO Auto-generated method stub

	}
	
	public static void data() throws Exception{
		//create a DNA alphabet
		DNAAlphabet dna = DNAAlphabet.SINGLETON;
		
		//define an arbitrary discrete alphabet
		DiscreteAlphabet discrete = new DiscreteAlphabet( false, "W", "S", "w", "x" );
		
		//define a discrete alphabet of numerical values
		DiscreteAlphabet numerical = new DiscreteAlphabet( 3, 10 );
		
		//define a continuous alphabet
		ContinuousAlphabet continuousInf = new ContinuousAlphabet();
		ContinuousAlphabet continuousPos = new ContinuousAlphabet( 0.0, 100.0 );
		
		//define a complementary alphabet
		GenericComplementableDiscreteAlphabet complementable = new GenericComplementableDiscreteAlphabet( true, new String[]{"A","B"}, new int[]{1,0} );
		
		//create an alphabet container for DNA
		AlphabetContainer dnaContainer = DNAAlphabetContainer.SINGLETON;
		
		//create a continuous alphabet container
		AlphabetContainer contContainer = new AlphabetContainer( continuousInf );
		
		//create a mixed alphabet container
		AlphabetContainer mixedContainer = new AlphabetContainer(dna, discrete, continuousPos);
		
		//create a section-wise defined alphabet container
		AlphabetContainer complex = new AlphabetContainer( new Alphabet[]{dna,continuousInf}, new int[]{0,0,1,0,1,1} );		
		
		//create a DNA-sequence from a string
		Sequence dnaSeq = Sequence.create( dnaContainer, "ACGTACGTACGT" );
		
		//create a continuous sequence from a string
		Sequence contSeq = Sequence.create( contContainer, "0.5 1.32642 99.5 20.4 5 7.7" , " " );
		
		//create a mixed sequence from a string
		Sequence mixedSeq = Sequence.create( mixedContainer, "C;x;5.67" , ";" );
		
		//create a sparse sequence
		Sequence sparse = new SparseSequence( dnaContainer, "ACGTACGTACGT" );
		
		//get length of sequence
		int length = dnaSeq.getLength();
		
		//obtain a single value
		int value = dnaSeq.discreteVal( 2 );
		double value2 = contSeq.continuousVal( 5 );
		
		//obtain a sub-sequence
		Sequence contSub = contSeq.getSubSequence( 2, 3 );
		
		//obtain a reverse complementary sequence
		Sequence revComp = dnaSeq.reverseComplement();
		
		//obtain the complement of a sub-sequence
		Sequence subComp = dnaSeq.complement( 3, 6 );
		
		//create a permuted sequence
		PermutedSequence permuted = new PermutedSequence( dnaSeq );
		
		//add annotation to a sequence
		Sequence annotatedDnaSeq = dnaSeq.annotate( true, new MotifAnnotation( "new motif", 3, 5, Strand.FORWARD ) );
		
		//retrieve annotation of a sequence
		SequenceAnnotation[] allAnnotations = annotatedDnaSeq.getAnnotation();
		
		//retrieve a specific annotation
		MotifAnnotation motif = (MotifAnnotation) annotatedDnaSeq.getSequenceAnnotationByType( "Motif", 0 );
		
		//create a data set from a DNA-FastA
		DNADataSet dnaDataSet = new DNADataSet( "myfile.fa" );
		
		//create a data set from a continuous tab-separated file
		DataSet contDataSet = new DataSet( contContainer, new SparseStringExtractor( "myfile.tab", '#' ), "\t" );
		
		//create a data set of sparse sequences
		DataSet sparseDataSet = SparseSequence.getDataSet( dnaContainer, new SparseStringExtractor( "myfile.fa", '>' ) );
		
		//retrieve a sequence from a sample
		Sequence fifth = dnaDataSet.getElementAt( 5 );
		
		for(int i=0;i<dnaDataSet.getNumberOfElements();i++){
			System.out.println(dnaDataSet.getElementAt( i ));
		}
		
		//iterate over sample
		for(Sequence seq : contDataSet){
			System.out.println(seq.getLength());
		}
		
		//infix sample
		DataSet infix = dnaDataSet.getInfixDataSet( 3, 10 );
		
		//suffix sample
		DataSet suffix = dnaDataSet.getSuffixDataSet( 7 );
		
		//complementary sample
		DataSet allRevComplements = dnaDataSet.getReverseComplementaryDataSet();
		
		//partition
		DataSet[] fiveParts = dnaDataSet.partition( 5, PartitionMethod.PARTITION_BY_NUMBER_OF_ELEMENTS );
		DataSet[] randParts = dnaDataSet.partition( PartitionMethod.PARTITION_BY_NUMBER_OF_SYMBOLS, 0.1, 0.2, 0.7 );
		
		//sub-sequences of user-defined length
		DataSet sliding = new DataSet( dnaDataSet, 8 );
		
		//create a data set from DNA-FastA and parse comment line
		DNADataSet dnaWithComments = new DNADataSet( "myfile.fa", '>', new SimpleSequenceAnnotationParser() );
		
		//obtain comment line from sequence
		String comment = dnaWithComments.getElementAt( 0 ).getAnnotation()[0].getResultAt( 0 ).getValue().toString();
		
		//create data set from DNA-FastA and parse entries of comment line
		DNADataSet dnaWithParsedComments = new DNADataSet( "myfile.fa", '>', new SplitSequenceAnnotationParser("=",";") );
		
		//obtain entries as annotation
		SequenceAnnotation[] allAnnotations2 = dnaWithParsedComments.getElementAt( 0 ).getAnnotation();
		
		//use BioJava to obtain sequences from genbank
		GenbankRichSequenceDB db = new GenbankRichSequenceDB();
		HashSet<String> idSet = new HashSet<String>( 2 );
		idSet.add( "NC_001284.2" );
		idSet.add( "NC_000932.1" );
		RichSequenceDB subDB = db.getRichSequences( idSet );
		RichSequenceIterator dbIterator = subDB.getRichSequenceIterator();
		
		DataSet fromBioJava = BioJavaAdapter.sequenceIteratorToDataSet( dbIterator, null );
		
		SequenceIterator backFromJstacs = BioJavaAdapter.dataSetToSequenceIterator( fromBioJava, true );
		
	}
	
	public static void xmlParser() throws Exception {
		
		StringBuffer buffer = new StringBuffer();
		
		//create and store primities
		int integer = 5;
		XMLParser.appendObjectWithTags( buffer, integer, "integer" );
		String bar = "hello world";
		XMLParser.appendObjectWithTags( buffer, bar, "foo" );
		
		//create and store arrays of primitives
		double[][] da = new double[4][6];
		XMLParser.appendObjectWithTags( buffer, da, "da" );
		
		//create and store Storable
		HomogeneousMM hMM = new HomogeneousMM( new HomMMParameterSet( DNAAlphabetContainer.SINGLETON, 4, "hmm(0)", (byte) 0 ) );
		XMLParser.appendObjectWithTags( buffer, hMM, "hMM" );
		
		//create and store arrays of Storables
		Storable[] storAr = ArrayHandler.createArrayOf( hMM, 5 );
		XMLParser.appendObjectWithTags( buffer, storAr, "storAr" );
		
		//parse primitives
		integer = (Integer) XMLParser.extractObjectForTags( buffer, "integer" );
		
		//parse arrays of primitives
		da = XMLParser.extractObjectForTags( buffer, "da", double[][].class );
		
		//parse Storables
		hMM = XMLParser.extractObjectForTags( buffer, "hMM", HomogeneousMM.class );
		
		//parse arrays of Storables
		storAr = (Storable[]) XMLParser.extractObjectForTags( buffer, "storAr" );
		
		//load specific type
		HomogeneousMM[] hmAr = ArrayHandler.createArrayOf( hMM, 5 );
		XMLParser.appendObjectWithTags( buffer, hmAr, "hmAr" );
		hmAr = (HomogeneousMM[]) XMLParser.extractObjectForTags( buffer, "hmAr" );
		
	}
	
	public static void parameters() throws Exception{
		
		//create simple parameter
		SimpleParameter simplePar = new SimpleParameter( DataType.INT, "Sequence length", "The required length of a sequence", true, new NumberValidator<Integer>( 1, 100 ), 10 );
		SimpleParameter simplePar2 = new SimpleParameter( DataType.STRING, "Name", "The name of the game", false );
		
		//create collection parameter from enum
		EnumParameter enumpar = new EnumParameter( DataType.class, "Data types", true );
		
		//create collection parameter
		SelectionParameter collPar = new SelectionParameter( DataType.DOUBLE, new String[]{"small", "large"}, new Double[]{5.0,5E6}, "Numbers", "A selection of numbers", true );
		
		//create collection parameter from subclasses
		collPar = SubclassFinder.getSelectionParameter( SequenceScoringParameterSet.class, "de", "Sequence scores", "All Sequence scores in Jstacs that can be created from parameter sets", true );
		
		//create simple parameter set
		SimpleParameterSet parSet = new SimpleParameterSet( simplePar,collPar );
		
		//implement own parameter set
		InstanceParameterSet myInstancePS = new InstanceParameterSet(TrainableStatisticalModel.class) {
			
			{
				this.initParameterList();
				this.parameters.add( new EnumParameter( DataType.class, "Data types", true ) );
				this.parameters.add( new SimpleParameter( DataType.INT, "Sequence length", "The accepted length of a sequence", true, new NumberValidator<Integer>( 1, 100 ), 10 )  );
			}
			
			@Override
			public String getInstanceName() {
				return "MyModel";
			}
			
			@Override
			public String getInstanceComment() {
				return "A fancy model created by me";
			}
		};
		
		//parameter set container
		ParameterSetContainer container = new ParameterSetContainer( "Set", "A set of parameters", parSet );
	}
	
	public static void results() throws Exception {
		
		//create numerical result
		NumericalResult res = new NumericalResult( "A double result", "This result contains some double value", 5.0 );
		
		//create categorical result
		CategoricalResult catRes = new CategoricalResult( "A boolean result", "This result contains some boolean", true );
		
		//create result set
		ResultSet resSet = new ResultSet( new Result[]{res,catRes} );
		
		//create and fill mean result set
		MeanResultSet mrs = new MeanResultSet();
		
		Random r = new Random();
		for(int i=0;i<10;i++){
			mrs.addResults( new NumericalResultSet( new NumericalResult( "Single", "A single result to be aggregated", r.nextDouble() ) ) );
		}
		System.out.println( mrs.getStatistics() );
		
	}
	
	public static void trainSMs() throws Exception {
		
		AlphabetContainer alphabet = DNAAlphabetContainer.SINGLETON;
		DataSet ds = new DNADataSet( "myfile.fa" );		
		
		//create models using model factory
		TrainableStatisticalModel pwm = TrainableStatisticalModelFactory.createPWM( alphabet, 10, 4.0 );
		TrainableStatisticalModel imm = TrainableStatisticalModelFactory.createInhomogeneousMarkovModel( alphabet, 12, 4.0, (byte) 2 );
		TrainableStatisticalModel pmm = TrainableStatisticalModelFactory.createPermutedMarkovModel( alphabet, 7, 4.0, (byte) 1 );
		TrainableStatisticalModel hmm = TrainableStatisticalModelFactory.createHomogeneousMarkovModel( alphabet, 400.0, (byte) 3 );
		TrainableStatisticalModel zoops = TrainableStatisticalModelFactory.createZOOPS( pwm, hmm, new double[]{4,4}, false );
		
		//train the PWM
		pwm.train( ds );
		
		//create hMM directly
		HomogeneousMM hmm2 = new HomogeneousMM( new HomMMParameterSet( alphabet, 4.0, "hmm(0)", (byte) 0 ) );
		
		//create BNM directly
		BayesianNetworkTrainSM bnm = new BayesianNetworkTrainSM( new BayesianNetworkTrainSMParameterSet( alphabet, 8, 4.0, "Bayesian network", ModelType.BN, (byte) 1, LearningType.ML_OR_MAP ) );
		
		//use hmm factory
		HMMTrainingParameterSet trainingPars = new BaumWelchParameterSet( 5, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 ),2 );
		Emission[] emissions = new Emission[]{new DiscreteEmission( alphabet, 4.0 ),new DiscreteEmission( alphabet, new double[]{2.0,1.0,1.0,2.0} )};
		AbstractHMM myHMM = HMMFactory.createErgodicHMM( trainingPars, 1, 4.0, 0.1, 100.0, emissions );
		
		//create hmm directly
		HigherOrderHMM hohmm = new HigherOrderHMM( trainingPars, new String[]{"A","B"}, emissions, 
				new TransitionElement( null, new int[]{0}, new double[]{4.0} ), 
				new TransitionElement( new int[]{0}, new int[]{0,1}, new double[]{2.0,2.0} ),
				new TransitionElement( new int[]{1}, new int[]{0}, new double[]{4.0} ));
		
		//print graphviz representation of HMM
		System.out.println( hohmm.getGraphvizRepresentation( null ) );
		
		//create mixture model of two PWMs, learn by EM
		MixtureTrainSM mixEm = new MixtureTrainSM( 8, new TrainableStatisticalModel[]{pwm,pwm}, 3, new double[]{4,0,4.0}, 1, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 ), Parameterization.LAMBDA );
		
		//create mixture model of two PWMs, learn by Gibbs-sampling
		MixtureTrainSM mixGibbs = new MixtureTrainSM( 8, new TrainableStatisticalModel[]{pwm,pwm}, 3, new double[]{4,0,4.0}, 100, 1000, new VarianceRatioBurnInTest( new VarianceRatioBurnInTestParameterSet( 3, 1.2 ) ) );
		
		//create strand model of iMM
		StrandTrainSM strandModel = new StrandTrainSM( imm, 3, 0.5, 1, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 ), Parameterization.LAMBDA );
		
		//create ZOOPS, learn by EM
		ZOOPSTrainSM zoops2 = new ZOOPSTrainSM( pwm, hmm, true, 4, 0.7, null, 1, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 ), Parameterization.LAMBDA );
		
		//implement own model using AbstractModel
		//TODO
	}
	
	public static void scoringFunctions() throws Exception {
		
		AlphabetContainer alphabet = DNAAlphabetContainer.SINGLETON;
		DataSet[] data = null;//TODO create
		
		//create BNSF
		BayesianNetworkDiffSM bnDsm = new BayesianNetworkDiffSM( new BayesianNetworkDiffSMParameterSet( alphabet, 10, 4.0, true, new BTExplainingAwayResidual( new double[]{4.0,4.0} ) ) );
		
		//create MMSF
		MarkovModelDiffSM mmDsm = new MarkovModelDiffSM( alphabet, 8, 4.0, true, new InhomogeneousMarkov( 1 ) );
		
		//train mmsfs discriminatively
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet( alphabet, 8, Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1, false, KindOfParameter.PLUGIN, true, 4 );
		GenDisMixClassifier cl = new GenDisMixClassifier( params, new CompositeLogPrior(), LearningPrinciple.MSP, mmDsm, mmDsm );
		cl.train( data );
		System.out.println(cl);
		
		//create HMMSF
		HomogeneousMMDiffSM hmmDsm = new HomogeneousMMDiffSM( alphabet, 3, 4.0, 100 );
		
		//create DHOMM
		DifferentiableEmission[] emissions = new DifferentiableEmission[]{new DiscreteEmission( alphabet, 4.0 ),new DiscreteEmission( alphabet, new double[]{2.0,1.0,1.0,2.0} )};
		NumericalHMMTrainingParameterSet trainingParameterSet = new NumericalHMMTrainingParameterSet( 3, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 ), 2, Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1 );
		DifferentiableHigherOrderHMM hmm = new DifferentiableHigherOrderHMM( trainingParameterSet, new String[]{"A","B"} , new int[]{0,1}, new boolean[]{true,true},emissions, true,4.0,
				new TransitionElement( null, new int[]{0}, new double[]{4.0} ), 
				new TransitionElement( new int[]{0}, new int[]{0,1}, new double[]{2.0,2.0} ),
				new TransitionElement( new int[]{1}, new int[]{0}, new double[]{4.0} ));
		
		//create mixture scoring function
		MixtureDiffSM mixDsm = new MixtureDiffSM( 3, true, mmDsm,mmDsm );
		
		//create strand scoring function
		
		StrandDiffSM strandDsm = new StrandDiffSM( mmDsm, 0.5, 1, true, InitMethod.INIT_BOTH_STRANDS );
		
		//create HiddenMotifsMixture
		ExtendedZOOPSDiffSM zoops = new ExtendedZOOPSDiffSM( ExtendedZOOPSDiffSM.CONTAINS_SOMETIMES_A_MOTIF, 500, 4, false, hmmDsm, strandDsm, null, true );
		
		//use IPSF
		
		IndependentProductDiffSM ipsf = new IndependentProductDiffSM( 4.0, true, bnDsm,mmDsm );
		
		//use NormalizedScoringFunctionModel
		
		DifferentiableStatisticalModelWrapperTrainSM trainSm = new DifferentiableStatisticalModelWrapperTrainSM( mmDsm, 4, Optimizer.QUASI_NEWTON_BFGS, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 ), 1E-6, 1 );
		trainSm.train( data[0] );
		
		//implement AbstractNormalizableScoringFunction
		//TODO
		
	}
	
	
	public static void optimization() throws Exception{
		
		//create numerical differentiable function
		NumericalDifferentiableFunction ndf = new NumericalDifferentiableFunction(1E-10) {
			
			@Override
			public int getDimensionOfScope() {
				return 2;
			}
			
			@Override
			public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException {
				return x[0]*x[0] + x[1]*x[1];
			}
		};
		
		//create differentiable function
		DifferentiableFunction df = new DifferentiableFunction() {
			
			@Override
			public int getDimensionOfScope() {
				return 2;
			}
			
			@Override
			public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException {
				return x[0]*x[0] + x[1]*x[1];
			}
			
			@Override
			public double[] evaluateGradientOfFunction( double[] x ) throws DimensionException, EvaluationException {
				return new double[]{2.0*x[0], 2.0*x[1]};
			}
		};
		
		
		//create termination condition
		AbstractTerminationCondition tc = new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 );
		AbstractTerminationCondition tc2 = new IterationCondition(100);
		
		//create combined termination condition
		TerminationCondition combined = new CombinedCondition( 2, tc, tc2 );
		
		//use optimizer
		double[] parameters = new double[df.getDimensionOfScope()];
		Optimizer.optimize( Optimizer.QUASI_NEWTON_BFGS, df, parameters, combined, 1E-6, new ConstantStartDistance( 1E-4 ), System.out );
		
	}
	
//classifier-section
	public static void classifier() throws Exception {		
		AlphabetContainer alphabet = DNAAlphabetContainer.SINGLETON;
		DataSet[] data = null;//TODO create
		
		//create models
		TrainableStatisticalModel pwm = TrainableStatisticalModelFactory.createPWM( alphabet, 10, 4.0 );
		
		//create and train model based classifier
		AbstractClassifier cl = new TrainSMBasedClassifier( pwm, pwm );
		
		//create the parameters for GenDisMixClassifier
		GenDisMixClassifierParameterSet ps = new GenDisMixClassifierParameterSet( alphabet, 10, (byte) 10, 1E-6, 1E-9, 1, false, KindOfParameter.PLUGIN, true, 2 );
		
		//create scoring functions
		DifferentiableStatisticalModel pwm2 = new BayesianNetworkDiffSM( alphabet, 10, 4.0, true, new InhomogeneousMarkov(0) );
		
		//create a GenDisMixClassifier using ML
		cl = new GenDisMixClassifier(ps, DoesNothingLogPrior.defaultInstance, LearningPrinciple.ML, pwm2, pwm2 );
		
		//create log prior
		LogPrior prior = new CompositeLogPrior();
		
		//create a GenDisMixClassifier using MSP
		cl = new GenDisMixClassifier(ps, prior, LearningPrinciple.MSP, pwm2, pwm2 );
		
		//create a GenDisMixClassifier using some hybrid learning principle
		cl = new GenDisMixClassifier(ps, prior, new double[]{0.4,0.1,0.5}, pwm2, pwm2 );
		
		//train
		cl.train( data );
		
		//classify sequence
		System.out.println( cl.classify( data[0].getElementAt(0) ) );
		
		//define performance measures
		PerformanceMeasureParameterSet measures = PerformanceMeasureParameterSet.createFilledParameters( false, 0.999, 0.95, 0.95, 1 );
		AbstractPerformanceMeasure[] m = {new AucROC(), new AucPR()};
		measures = new PerformanceMeasureParameterSet( m );
		
		//assess model based classifier on test data
		System.out.println( cl.evaluate( measures, true, data ) );
		
		//assess classifiers in CV
		NumericalPerformanceMeasureParameterSet numMeasures = PerformanceMeasureParameterSet.createFilledParameters();
		
		ClassifierAssessment assessment = new KFoldCrossValidation( cl );
		KFoldCrossValidationAssessParameterSet params = new KFoldCrossValidationAssessParameterSet( PartitionMethod.PARTITION_BY_NUMBER_OF_ELEMENTS, cl.getLength(), true, 10 );
		System.out.println( assessment.assess( numMeasures, params, data ) );				
	}
	
	public static void alignment(){
		
		//create costs
		Costs costs = new SimpleCosts( 0, 1, 1, 0.5 );
		
		Sequence seq1=null, seq2=null;
		//create alignment of two string
		Alignment align = new Alignment( AlignmentType.GLOBAL, costs );
		System.out.println( align.getAlignment( seq1, seq2 ) );
		
	}
	
	public void utils() throws Exception{
		
		//create REnvironment
		REnvironment re = new REnvironment();
		
		double[] values = new double[10];
		//use REnvironment
		re.createVector( "values", values );
		
		//plot
		re.plotToPDF( "plot(values,t=\"l\");", "values.pdf", true );
		
		double[][] twodim = new double[5][5];
		//ArrayHandler clone
		ArrayHandler.clone( twodim );
		
		//ArrayHandler create
		TrainableStatisticalModel pwm = TrainableStatisticalModelFactory.createPWM( DNAAlphabetContainer.SINGLETON, 10, 4.0 );
		TrainableStatisticalModel[] models = ArrayHandler.createArrayOf( pwm, 10 );
		
		//ArrayHandler cast
		BayesianNetworkTrainSM[] bns = ArrayHandler.cast( BayesianNetworkTrainSM.class, models );
		
		//ToBo.max
		double max = ToolBox.max( values );
		
		//ToBo.sum
		double sum = ToolBox.sum( values );
		
		//ToBo.maxIndex
		int maxIndex = ToolBox.getMaxIndex( values );
		
		//SafeOutputStream
		OutputStream stream = SafeOutputStream.getSafeOutputStream( System.out );
		
		//IntList
		
		//DoubleList
		
		//SubclassFinder
		LinkedList<Class<? extends TrainableStatisticalModel>> list = SubclassFinder.findInstantiableSubclasses( TrainableStatisticalModel.class, "de.jstacs" );
		
		//UserTime
		
	}
	

}
