package projects.dream2016.mix;
import java.util.Arrays;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NotTrainedException;
import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.DifferentiableFunction;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.NegativeDifferentiableFunction;
import de.jstacs.algorithms.optimization.NumericalDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.StartDistanceForecaster;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.MixtureDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.DifferentiableStatisticalModelWrapperTrainSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;

//TODO private!?

/**
 * 
 * @author Jens Keilwagen
 */
public class NewMixtureClassifier extends AbstractScoreBasedClassifier implements OptimizableClassifier {

	private AbstractScoreBasedClassifier[] component;
	private MyMixtureScoringFunction model;
	private LogPrior mixPrior;
	private Vote vote;
	private Training training;
	private int threads;
	private int starts;
	
	private OptimizableClassifier[] optComponent;
	private double[] helpArray, logClassifierProbs, logMixScores, params, mixParams, mixGrad;
	private IntList[] indices;
	private DoubleList[] partialDer;
		
	private byte algo = (byte) 10;//Optimizer.QUASI_NEWTON_BFGS;
	private double eps = 1E-7;//TODO
	private double linEps = 1E-7;
	
	public NewMixtureClassifier( int threads, Training training, int starts, DifferentiableStatisticalModel[] componentSF,
			AbstractScoreBasedClassifier[] componentClassifiers, Vote vote, LogPrior prior ) throws Exception {
		super(componentClassifiers[0].getAlphabetContainer(), componentClassifiers[0].getLength(), componentClassifiers[0].getNumberOfClasses() );
		if( componentSF != null && componentSF.length != componentClassifiers.length ) {
			throw new IllegalArgumentException();
		}
		helpArray = new double[getNumberOfClasses()];
		model = new MyMixtureScoringFunction(training==Training.COMBINED?1:starts,true,ArrayHandler.clone(componentSF));
		if( starts <= 0 ) {
			throw new IllegalArgumentException("The number of starts has to be positive.");
		} else {
			this.starts = starts;
		}		
		component = ArrayHandler.clone(componentClassifiers);
		logClassifierProbs = new double[componentClassifiers.length];
		logMixScores = new double[componentSF.length];
		this.training = training;
		if( training == Training.COMBINED ) {
			optComponent = ArrayHandler.cast( OptimizableClassifier.class, component );
		} else {
			optComponent = null;
		}
		this.threads = threads;
		mixPrior = prior;
		this.vote = vote;
	}
	
	
	public NewMixtureClassifier( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#extractFurtherClassifierInfosFromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		super.extractFurtherClassifierInfosFromXML( xml );
		model = XMLParser.extractObjectForTags( xml, "model", MyMixtureScoringFunction.class );
		component = XMLParser.extractObjectForTags( xml, "componentClassifier", AbstractScoreBasedClassifier[].class );
		StringBuffer pr = XMLParser.extractForTag( xml, "mixPrior" );
		if( pr != null )
		{
			String className = XMLParser.extractObjectForTags( pr, "className", String.class );
			try
			{
				mixPrior = (LogPrior) Class.forName( className ).getConstructor( new Class[]{ StringBuffer.class } )
						.newInstance( pr );
			}
			catch( NoSuchMethodException e )
			{
				NonParsableException n = new NonParsableException( "You must provide a constructor " + className
						+ "(StringBuffer)." );
				n.setStackTrace( e.getStackTrace() );
				throw n;
			}
			catch( Exception e )
			{
				NonParsableException n = new NonParsableException( "problem at " + className + ": " + e.getMessage() );
				n.setStackTrace( e.getStackTrace() );
				throw n;
			}
		}
		else
		{
			mixPrior = DoesNothingLogPrior.defaultInstance;
		}
		vote = XMLParser.extractObjectForTags( xml, "vote", Vote.class );
		training = XMLParser.extractObjectForTags( xml, "training", Training.class );
		threads = XMLParser.extractObjectForTags( xml, "threads", int.class );
		starts = XMLParser.extractObjectForTags( xml, "starts", int.class );
		if( training == Training.COMBINED ) {
			optComponent = ArrayHandler.cast( OptimizableClassifier.class, component );
			try {
				reset();
			} catch( Exception e ) {
				throw new NonParsableException( e.getMessage() );
			}
		} else {
			optComponent = null;
		}
		helpArray = new double[getNumberOfClasses()];
		logClassifierProbs = new double[component.length];
		logMixScores = new double[logClassifierProbs.length];
		
	}
	
	/* (non-Javadoc)
	 * @see de.jstacs.classifier.AbstractClassifier#getFurtherClassifierInfos()
	 */
	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = super.getFurtherClassifierInfos();
		XMLParser.appendObjectWithTags( xml, model, "model" );
		XMLParser.appendObjectWithTags( xml, component, "componentClassifier" );
		if( !(mixPrior instanceof DoesNothingLogPrior) )
		{
			StringBuffer pr = new StringBuffer( 1000 );
			pr.append( "<mixPrior>\n" );
			XMLParser.appendObjectWithTags( pr, mixPrior.getClass().getName(), "className" );
			pr.append( mixPrior.toXML() );
			pr.append( "\t</mixPrior>\n" );
			xml.append( pr );
		}
		XMLParser.appendObjectWithTags( xml, vote, "vote" );
		XMLParser.appendObjectWithTags( xml, training, "training" );
		XMLParser.appendObjectWithTags( xml, threads, "threads" );
		XMLParser.appendObjectWithTags( xml, starts, "starts" );
		return xml;
	}
	
	public NewMixtureClassifier clone() throws CloneNotSupportedException {
		NewMixtureClassifier clone = (NewMixtureClassifier) super.clone();
		clone.component = ArrayHandler.clone( component );
		if( optComponent != null ) {
			clone.optComponent = ArrayHandler.cast( OptimizableClassifier.class, clone.component );
		}
		clone.model = (MyMixtureScoringFunction) model.clone();
		clone.mixPrior = mixPrior.getNewInstance();
		clone.logClassifierProbs = logClassifierProbs.clone();
		clone.logMixScores = logMixScores.clone();
		if( params != null ) {
			clone.params = params.clone();
			clone.mixParams = mixParams.clone();
			clone.mixGrad = mixGrad.clone();
			clone.indices = new IntList[indices.length];
			clone.partialDer = new DoubleList[partialDer.length];
			for( int i = 0; i < optComponent.length; i++ ) {
				clone.indices[i] = new IntList();
				clone.partialDer[i] = new DoubleList();
			}
		}
		clone.helpArray = helpArray.clone();
		return clone;
	}
	
	public void setVote( Vote v ) {
		vote = v;
	}
	
	@Override
	protected double getScore(Sequence seq, int i, boolean check)
			throws IllegalArgumentException, NotTrainedException, Exception {
		return getLogProb( seq, i, check, vote );
	}
	
	private double getLogProb(Sequence seq, int i, boolean check, Vote vote)
		throws EvaluationException {
		try {
			if( check ) {
				check( seq );
			}
		
			int k, j;
			switch (vote) {
			case DOC:
				k = model.getIndexOfMaximalComponentFor(seq, 0); 
				for( j = 0; j < helpArray.length; j++ ) {
					helpArray[j] = component[k].getScore(seq, j);
				}
				return helpArray[i] - Normalisation.getLogSum(helpArray);
			case VOC:
				model.getLogScores( seq, logMixScores ); 
				double logSum = Normalisation.getLogSum( logMixScores );
				for( k = 0; k < logMixScores.length; k++ ) {
					for( j = 0; j < helpArray.length; j++ ) {
						helpArray[j] = component[k].getScore(seq, j);
					}
					logMixScores[k] += helpArray[i] - Normalisation.getLogSum(helpArray);
				}
				return Normalisation.getLogSum(logMixScores) - logSum;
			default:
				throw new RuntimeException();
			}
		} catch( Exception e ) {
			EvaluationException ee = new EvaluationException( e.getClass() + ": "+  e.getMessage() );
			ee.setStackTrace( e.getStackTrace() );
			throw ee;
		}
	}

	@Override
	public CategoricalResult[] getClassifierAnnotation() {
		//TODO
		CategoricalResult[] res = new CategoricalResult[getNumberOfClasses() + 1];
		res[0] = new CategoricalResult( "classifier", "a <b>short</b> description of the classifier", getInstanceName() );
		for( int i = 1; i < res.length; i++ ) {
			res[i] = new CategoricalResult( "class info " + (i-1), "some information about the class", "" );//TODO
		}
		return res;
	}

	@Override
	public String getInstanceName() {
		StringBuffer sb = new StringBuffer();
		sb.append( getClass().getSimpleName() + "(model: " + model.getInstanceName() + "; classifier: " );
		sb.append( component[0].getInstanceName() );
		for( int i = 1; i < component.length; i++ ) {
			sb.append( ", " + component[i].getInstanceName() );
		}
		sb.append( ")" );
		return sb.toString();
	}

	@Override
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		LinkedList<NumericalResult> list = new LinkedList<NumericalResult>();
		NumericalResultSet nrs;
		/*TODO name clash
		for( int i = 0; i < component.length; i++ ) {
			nrs = component[i].getNumericalCharacteristics();
			if( nrs.getNumberOfResults() > 0 ) {
				list.add( new NumericalResult("classifier " + i,"the index of the classifier",i) );
				for( int j = 0;j < nrs.getNumberOfResults(); j++ ) {
					list.add( nrs.getResultAt( j ) );
				}
			}
		}
		*/
		return new NumericalResultSet( list );
	}

	private static final String XML_TAG = "MixtureClassifier";
	
	@Override
	protected String getXMLTag() {
		return XML_TAG;
	}

	@Override
	public boolean isInitialized() {
		int i = 0;
		while( i < component.length && component[i].isInitialized() ) {
			i++;
		}
		return i == component.length && model.isInitialized();
	}

	
	private DataHandler init( DataSet[] s, double[][] weights, boolean trainModel ) throws Exception {
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		DoubleList w = new DoubleList();
		for( int j = 0; j < s.length; j++ ) {
			for( int n = 0; n < s[j].getNumberOfElements(); n++ ) {
				seqs.add(s[j].getElementAt(n) );
				if( weights != null && weights[j] != null ) {
					w.add( weights[j][n] );
				}
			}
		}

		DataSet[] all = { new DataSet( "all", seqs.toArray( new Sequence[0] ) ) };
		double[][] weight = { (w.length() == 0 ? null : w.toArray()) };
		if( model.getNumberOfParameters()>model.getNumberOfComponents() && trainModel ) {
			DifferentiableStatisticalModelWrapperTrainSM myModel = new DifferentiableStatisticalModelWrapperTrainSM(model,threads,algo,new SmallDifferenceOfFunctionEvaluationsCondition(1E-3/*TODO*/),linEps,starts, this.mixPrior.getNewInstance());
			//XXX myModel.setOutputStream(null);
			myModel.train( all[0], weight[0] );
			model = (MyMixtureScoringFunction) myModel.getFunction();
		} else {
			model.initializeFunction( 0, false, all, weight );
		}
/*
		double[] p = model.getCurrentParameterValues();
		Arrays.fill(p, 0,model.getNumberOfComponents(), 0);
		p[0]=4;
		model.setParameters(p, 0);/**/
System.out.println(model);
/**/
		weight = null;
		all = null;
		
		DataHandler dh = splitData(model, s, weights, training==Training.SEPARATELY_DOC );
		
		if( optComponent != null ) {
			for( int i = 0; i < optComponent.length; i++ ) {
				optComponent[i].initialize( dh.getSplits(i), dh.getWeightsOfSplits(i)  );
			}
		}
		return dh;
	}
	
	@Override
	public void train(DataSet[] s, double[][] weights) throws Exception {
		if( training == Training.COMBINED ) {
			if( weights == null ){
				weights = new double[s.length][];
			}
			for( int i = 0; i < s.length; i++ ) {
				if( weights[i] == null ) {
					weights[i]= new double[s[i].getNumberOfElements()];
					Arrays.fill( weights[i], 1 );
				}
			}
			
			NewMixtureClassifier bestClone = this;
			MSPClassifierObjective obj;
			DifferentiableFunction neg;
			double[] params;
			double current, best = Double.NEGATIVE_INFINITY;
			
			DataSet[] sampled = new DataSet[s.length];
			double[][] sampledWeights = new double[s.length][];
			Pair<DataSet[], double[][]> p;
			TerminationCondition term = new SmallDifferenceOfFunctionEvaluationsCondition( eps );
			StartDistanceForecaster sd = new ConstantStartDistance(1);
			boolean trained = false;
			Exception last = null;
			for( int i = 0; i < starts; i++ ) {
System.out.println(i + " ========================");
				try {
					this.initializeRandomly(); //this should clear all parameters
					
					for( int j = 0; j < sampled.length; j++ ) {
						p = s[j].partition( weights[j], PartitionMethod.PARTITION_BY_NUMBER_OF_ELEMENTS, 0.95, 0.05 );
						sampled[j] = p.getFirstElement()[1];
						sampledWeights[j] = p.getSecondElement()[1];
					}
					
					DataHandler dh = init(sampled, sampledWeights, true);
					dh = splitData( model, s, weights, false );
					for( int j = 0; j < component.length; j++ ) {
						component[j].train( dh.getSplits(j), dh.getWeightsOfSplits(j) );
					}

					//clone.initializeRandomly();
					
					//We need new MSPClassifierObjectives here, because in case of an Exception,
					//stopThreads is called, workers is set to null and subsequent optimizations
					//(and calls to stopThreads) must fail
					obj = new MSPClassifierObjective(threads,this,s,weights,true);
					neg = new NegativeDifferentiableFunction(obj);
					
					obj.reset();
					
					params = obj.getParameters( KindOfParameter.PLUGIN );
					/*
					obj.setDataAndWeights(sampled, sampledWeights);
					Optimizer.optimize(algo, neg, params, term, linEps, sd, System.out );
					obj.setDataAndWeights(s, weights);
	*/
					Optimizer.optimize(algo, neg, params, term, linEps, sd, System.out );
					current = obj.evaluateFunction( params );
	System.out.println( "start " + i + ":\t" + current );
					if( current > best ) {
						bestClone = this.clone();
						best = current;
					}
					trained = true;
					

					obj.stopThreads();
					
				} catch( Exception e ) {
					last = e;
	System.out.println( "An exception was thrown. " + e.getMessage() );
				}
			}
			
			if( !trained ) {
				throw last;
			}
			this.component = bestClone.component;
			bestClone.component = null;
			optComponent = ArrayHandler.cast( OptimizableClassifier.class, component );
			this.model = bestClone.model;
			bestClone.model = null;
			this.params = null;
			reset();
		} else {
			DataHandler dh = init( s, weights, true );
			for( int i = 0; i < component.length;i ++ ) {
				DataSet[] sp = dh.getSplits(i);
System.out.println( i + ")\t" + sp[0].getNumberOfElements() + " vs. " + sp[1].getNumberOfElements() );
				component[i].train( sp, dh.getWeightsOfSplits(i) );
			}
		}
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append( model );
		for( int i = 0; i < component.length; i++ ) {
			sb.append( "\nclassifier " + i + ":\n" );
			sb.append( component[i] );
		}
		return sb.toString();
	}
	
	public static DataHandler splitData( MyMixtureScoringFunction mix, DataSet[] data, double[][] weights, boolean doc ) throws EmptyDataSetException, WrongAlphabetException {
		int k = mix.getNumberOfComponents();
		LinkedList<Sequence>[][] seqList = new LinkedList[k][data.length];
		DoubleList[][] weightsList = new DoubleList[k][data.length];
		for( int j, i = 0; i < k; i++ ) {
			for( j = 0; j < data.length; j++ ) {
				seqList[i][j] = new LinkedList<Sequence>();
				weightsList[i][j] = new DoubleList();
			}
		}
		
		Sequence seq;
		double w;
		double[] scores = new double[mix.getNumberOfComponents()]; 
		for( int anz, n, d = 0; d < data.length; d++ ) {
			anz = data[d].getNumberOfElements();
			w = 1;
			for( n = 0; n < anz; n++ ) {
				seq = data[d].getElementAt(n);
				if( weights != null && weights[d] != null ) {
					w = weights[d][n];
				}
				if( doc ) {
					k = mix.getIndexOfMaximalComponentFor(seq, 0);
					seqList[k][d].add( seq );
					weightsList[k][d].add( w );
				} else {
					seqList[0][d].add( seq );
					mix.getLogScores( seq, scores );
					Normalisation.logSumNormalisation( scores );
					for( k = 0; k < scores.length; k++ ) {
						weightsList[k][d].add( w*scores[k] );
					}
				}
				
			}
		}
		if( doc ) {
			return new DataHandler( weightsList, seqList );
		} else {
			return new DataHandler( weightsList, seqList[0] );
		}
		
	}
	
	public void addGradient( double[] grad, int start ) throws EvaluationException {
		Arrays.fill(  mixGrad, 0 );
		mixPrior.addGradientFor( mixParams, mixGrad );
		for( int i = 0; i < mixGrad.length; i++, start++ ) {
			grad[start] += mixGrad[i];
		}
		for( int i = 0; i < optComponent.length; i++ ) {
			optComponent[i].addGradient( grad, start );
			start += optComponent[i].getNumberOfParameters();
		}
	}

	public double[] getCurrentParameterValues( KindOfParameter kind ) throws Exception {
		return getCurrentParameterValues( kind, getNumberOfParameters(), null );
	}
	
	private double[] getCurrentParameterValues( KindOfParameter kind, int numOfParams, double[] arrayToFill ) throws Exception {
		if( arrayToFill == null ) {
			arrayToFill = new double[params.length];
		}
		mixParams = model.getCurrentParameterValues();
		numOfParams = mixParams.length;
		System.arraycopy( mixParams, 0, arrayToFill, 0, numOfParams );
		double[] part;
		for( int i = 0; i < optComponent.length; i++) {
			part = optComponent[i].getCurrentParameterValues( kind );
			System.arraycopy( part, 0, arrayToFill, numOfParams, part.length );
			numOfParams += part.length;
		}
		return arrayToFill;
	}

	public double getLogPriorTerm() throws DimensionException, EvaluationException {
		double logPriorTerm = mixPrior.evaluateFunction( mixParams );
		for( int i = 0; i < optComponent.length; i++ ) {
			logPriorTerm += optComponent[i].getLogPriorTerm();
		}	
		return logPriorTerm;
	}

	public double getLogProb( int classIndex, Sequence seq ) throws EvaluationException {
		return getLogProb(classIndex, seq,false);
	}
	
	private double getLogProb( int classIndex, Sequence seq, boolean out ) throws EvaluationException {
		//return getLogProb( seq, classIndex, false, Vote.VOC );
		model.getLogScores( seq, logMixScores );
		double logSum = Normalisation.getLogSum( logMixScores );
		for( int i = 0; i < optComponent.length; i++ ) {
			logClassifierProbs[i] = optComponent[i].getLogProb( classIndex, seq ) + logMixScores[i] - logSum;
		}
		return Normalisation.getLogSum( logClassifierProbs );
	}

//TODO check Jan
	public double getLogProbAndPartialDerivations( int classIndex, Sequence seq, IntList indices, DoubleList partialDer ) {
		//compute mixture weights
		model.getLogScores( seq, logMixScores );
		double nlms = Normalisation.getLogSum( logMixScores );
		
		//get all partial derivations
		for( int i = 0; i < optComponent.length; i++ ) {
			this.indices[i].clear();
			this.partialDer[i].clear();
			
			logClassifierProbs[i] = optComponent[i].getLogProbAndPartialDerivations( classIndex, seq, this.indices[i], this.partialDer[i] ) +logMixScores[i]-nlms;
		}
		
		//compute final result
		double result = Normalisation.logSumNormalisation( logClassifierProbs );

		//individual component classifier derivations
		int offset = mixParams.length;		
		for( int j, i = 0; i < this.indices.length; i++ ) {
			for( j = 0; j < this.indices[i].length(); j++ ) {
				indices.add( offset + this.indices[i].get( j ) );
				partialDer.add( logClassifierProbs[i] * this.partialDer[i].get( j ) );
			}
			offset += optComponent[i].getNumberOfParameters();
		}
		
		//mixture derivations
		model.getPartialDerivations( seq, logClassifierProbs, indices, partialDer );
		
		int i = 0;
		while( i < indices.length() ) {
			if( Double.isNaN( partialDer.get(i) ) ) {
				System.out.println("partial derivation became: NaN");
				System.out.println(classIndex + "\t" + seq);
				System.out.println(Arrays.toString(this.params) );
				System.out.println(i + "\t" + partialDer.get(i));
				System.exit(1);
				
			}
			i++;
		}

		return result;
	}

	public int getNumberOfParameters() {
		return params == null ? DifferentiableStatisticalModel.UNKNOWN : params.length;
	}

	public void initialize( DataSet[] s, double[][] weights ) throws Exception {
		init( s, weights, true );
	}
	
	public void initializeRandomly() throws Exception {
		model.initializeFunctionRandomly( false );
		for( int j = 0; j < optComponent.length; j++ ) {
			optComponent[j].initializeRandomly();
		}
		//System.out.println( this );
	}

	public void reset() throws Exception {
		mixPrior.set( true, model );
		for( int i = 0; i < optComponent.length; i++ ) {
			optComponent[i].reset();
		}
		int n, numOfParams = model.getNumberOfParameters();
		if( numOfParams != DifferentiableStatisticalModel.UNKNOWN ) {
			for( int i = 0; i < optComponent.length; i++ ) {
				n = optComponent[i].getNumberOfParameters();
				if( n == DifferentiableStatisticalModel.UNKNOWN ) {
					numOfParams = n;
					break;
				} else {
					numOfParams += n;
				}
			}
		}
		if( numOfParams != DifferentiableStatisticalModel.UNKNOWN ) {
			if( params == null || params.length != numOfParams ) {
				params = new double[numOfParams];
				mixParams = new double[model.getNumberOfParameters()];
				mixGrad = new double[mixParams.length];
			}
			getCurrentParameterValues( KindOfParameter.PLUGIN, numOfParams, params );
		} else {
			params = null;
		}
		if( indices == null ) {
			indices = new IntList[optComponent.length];
			partialDer = new DoubleList[optComponent.length];
			for( int i = 0; i < optComponent.length; i++ ) {
				indices[i] = new IntList();
				partialDer[i] = new DoubleList();
			}
		}
	}

	public void setParameters( double[] params, int start ) throws Exception {
		System.arraycopy( params, start, this.mixParams, 0, this.mixParams.length );
		System.arraycopy( params, start, this.params, 0, this.params.length );
		model.setParameters( params, start );
		start += model.getNumberOfParameters();
		for( int i = 0; i < optComponent.length; i++ ) {
			optComponent[i].setParameters( params, start );
			start += optComponent[i].getNumberOfParameters();
		}
	}
	
	
	/**
	 * 
	 * @author Jens Keilwagen
	 */
	public static class DataHandler {
		
		private DataSet[][] splits;
		private double[][][] weights;
		
		public DataHandler( DoubleList[][] weightsList, LinkedList<Sequence>[]... seqList ) throws EmptyDataSetException, WrongAlphabetException {
			int noOfSplits = weightsList.length, noOfClasses = weightsList[0].length;;
			splits = new DataSet[noOfSplits][];
			weights = new double[noOfSplits][][];
			Sequence[] empty = new Sequence[0];
			DataSet[] help = null;
			if( seqList.length == 1 ){
				help = new DataSet[noOfClasses];
				for( int c = 0; c < noOfClasses ; c++ ) {
					help[c] = new DataSet( "", seqList[0][c].toArray( empty ) );
				}
			}
			for( int c, s = 0; s < noOfSplits; s++ ) {
				splits[s] = new DataSet[noOfClasses];
				weights[s] = new double[noOfClasses][];
				for( c = 0; c < noOfClasses ; c++ ) {
					if( seqList.length == 1 ){
						splits[s][c] = help[c];
					} else {
						splits[s][c] = new DataSet( "", seqList[s][c].toArray( new Sequence[0] ) );
					}
					if( weightsList[s][c].length() == 0 ) {
						weightsList[s][c] = null;
					} else {
						weights[s][c] = weightsList[s][c].toArray();
					}
				}
			}
			
		}
		
		public int getNumberOfSplits() {
			return splits.length;
		}
		
		public DataSet[] getSplits( int i ) {
			return splits[i];
		}
		
		public double[][] getWeightsOfSplits( int i ) {
			return weights[i];
		}
	}
	
	/**
	 * This extension of {@link MixtureScoringFunction} provides additional methods that are needed for the {@link NewMixtureClassifier}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class MyMixtureScoringFunction extends MixtureDiffSM {

		public MyMixtureScoringFunction(int starts, boolean plugIn, DifferentiableStatisticalModel[] component)
				throws CloneNotSupportedException {
			super(starts, plugIn, component);
		}
		
		public MyMixtureScoringFunction( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * Fills the array with {@latex.inline $\\log score(k, x | \\lambda)$} for sequence x.
		 * 
		 * @param seq the sequence x
		 * @param logScores the array
		 */
		public void getLogScores( Sequence seq, double[] logScores ) {
			fillComponentScores(seq, 0);
			System.arraycopy( componentScore, 0, logScores, 0, componentScore.length );
		}
		
		public void getPartialDerivations( Sequence seq, double[] gamma, IntList indices, DoubleList partialDer ) {
			for( int i = 0; i < function.length; i++ ) {
				iList[i].clear();
				dList[i].clear();
				componentScore[i] = logHiddenPotential[i] + function[i].getLogScoreAndPartialDerivation( seq, 0, iList[i], dList[i] );
			}
			Normalisation.logSumNormalisation(componentScore);

			int k = paramRef.length-2, l = paramRef[k+1] - paramRef[k];
			for( int i = 0, j; i < function.length; i++ ) {
				double delta_i = gamma[i] - componentScore[i];
				//component model derivations
				for( j = 0; j < iList[i].length(); j++ ) {
					indices.add( paramRef[i] + iList[i].get(j) );
					partialDer.add( dList[i].get(j) * delta_i );
				}
				//mixture derivations
				if( i < l ) {
					indices.add( paramRef[k] + i );
					partialDer.add( delta_i );
				}
			}
		}
		
		protected boolean determineIsNormalized(){
			return false;
		}
		

		public void init( DataSet[] s, double[][] weights ) throws Exception {
			double[] stat = new double[function.length];
			SmallDifferenceOfFunctionEvaluationsCondition eps = new SmallDifferenceOfFunctionEvaluationsCondition(1E-11);
			for( int i = 0; i < stat.length; i++ ) {
				DifferentiableStatisticalModelWrapperTrainSM myModel = new DifferentiableStatisticalModelWrapperTrainSM(function[i],2,(byte) 10,eps,1E-9,1);
				myModel.train( s[i], weights==null?null:weights[i] );
				function[i] = myModel.getFunction();
				if( weights == null || weights[i] == null ) {
					stat[i] = s[i].getNumberOfElements();
				} else {
					throw new OperationNotSupportedException( "yet" );
				}
			}
			computeHiddenParameter( stat, true ); 
		}
	}
	
	/**
	 * Enum for the kind of decision.
	 * 
	 * @author Jens Keilwagen
	 */
	public static enum Vote {
		DOC, VOC
	}
	
	/**
	 * Enum for the kind of training.
	 * 
	 * @author Jens Keilwagen
	 */
	public static enum Training {
		SEPARATELY_DOC, SEPARATELY_VOC, COMBINED
	}
	
	/**
	 * This extension of {@link GenDisMixClassifier} implements the methods of {@link OptimizableClassifier} that allow to use this classifier in the {@link NewMixtureClassifier}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class OptimizableMSPClassifier extends GenDisMixClassifier implements OptimizableClassifier {

		private double[] parameter, grad, helpArray;
		private IntList[] indi;
		private DoubleList[] partDer;
		
		public OptimizableMSPClassifier( GenDisMixClassifierParameterSet params, LogPrior logPrior, DifferentiableStatisticalModel... score ) throws CloneNotSupportedException {
			super(params, logPrior, 0, LearningPrinciple.getBeta(LearningPrinciple.MSP), score);
			init();
		}
		
		public OptimizableMSPClassifier( StringBuffer xml ) throws NonParsableException {
			super( xml );
			init();
		}
		
		public OptimizableMSPClassifier clone() throws CloneNotSupportedException {
			OptimizableMSPClassifier clone = (OptimizableMSPClassifier) super.clone();
			if( parameter != null ) {
				clone.parameter = parameter.clone();
				clone.grad = grad.clone();
			}
			//recreate indi, partDer, helpArray
			clone.init();
			return clone;
		}
		
		private void init() {
			helpArray = new double[getNumberOfClasses()];
			indi = new IntList[helpArray.length];
			partDer = new DoubleList[helpArray.length];
			for( int i = 0; i < helpArray.length; i++ ) {
				indi[i] = new IntList();
				partDer[i] = new DoubleList();
			}
		}

		public void addGradient( double[] grad, int start ) throws EvaluationException {
			Arrays.fill( this.grad, 0 );
			prior.addGradientFor( parameter, this.grad );
			for( int i = 0; i < this.grad.length; i++ ) {
				grad[start+i] += this.grad[i]; 
			}
		}

		public double getLogPriorTerm() throws DimensionException, EvaluationException {
			return prior.evaluateFunction( parameter );
		}

		public double getLogProb( int classIndex, Sequence seq ) {
			for( int i = 0; i < score.length; i++ ) {
				helpArray[i] = getClassWeight( i ) + score[i].getLogScoreFor( seq );
			}
			return helpArray[classIndex] - Normalisation.getLogSum( helpArray );
		}

		public double getLogProbAndPartialDerivations( int classIndex, Sequence seq, IntList indices, DoubleList partialDer ) {
//TODO check Jan
			for( int i = 0; i < score.length; i++ ) {
				indi[i].clear();
				partDer[i].clear();
				helpArray[i] = getClassWeight( i ) + score[i].getLogScoreAndPartialDerivation( seq, indi[i], partDer[i] );
			}
			double logProb = helpArray[classIndex] - Normalisation.logSumNormalisation( helpArray );
			for( int i = 0, offset = score.length; i < score.length; i++ ) {
				indices.add( i );
				partialDer.add( (i==classIndex?1:0) - helpArray[i] );
				for( int j = 0; j < indi[i].length(); j++ ) {
					indices.add( offset + indi[i].get( j ) );
					partialDer.add( partDer[i].get( j ) * ( (i==classIndex?1:0) - helpArray[i] ) );
				}
				offset += score[i].getNumberOfParameters();
			}
			return logProb;
		}

		public int getNumberOfParameters() {
			int i = 0, num = 0, a;
			while( i < score.length ) {
				a = score[i++].getNumberOfParameters();
				if( a == DifferentiableStatisticalModel.UNKNOWN ) {
					return DifferentiableStatisticalModel.UNKNOWN;
				}
				num += a;
			}
			return getNumberOfClasses() + num;
		}

		public void initialize( DataSet[] data, double[][] weights) throws Exception {
			for( int i = 0; i < score.length; i++ ) {
				score[i].initializeFunction( i, false, data, weights );
			}			
			setClassWeights( false, new double[getNumberOfClasses()] );
			fillParameters();
		}
		
		public void initializeRandomly() throws Exception {
			for( int i = 0; i < score.length; i++ ) {
				score[i].initializeFunctionRandomly( false );
			}			
			setClassWeights( false, new double[getNumberOfClasses()] );
			fillParameters();
		}
		
		public void train( DataSet[] d, double[][] weight ) throws Exception {
			super.train(d, weight);
			fillParameters();
		}

		public void setParameters( double[] params, int start ) throws Exception {
			setClassWeights( false, params, start );
			start += getNumberOfClasses();
			for( int i = 0; i < score.length; i++ ) {
				score[i].setParameters( params, start );
				start += score[i].getNumberOfParameters();
			}
			fillParameters();
		}

		private void fillParameters() throws Exception {
			double[] help;
			int start = getNumberOfParameters();
			if( parameter == null || parameter.length != start ) {
				parameter = new double[start];
				grad = new double[start];
			}
			start = getNumberOfClasses();
			for( int i = 0; i < start; i++ ) {
				parameter[i] = getClassWeight( i );
			}
			for( int i = 0; i < score.length; i++ ) {
				help = score[i].getCurrentParameterValues();
				System.arraycopy( help, 0, parameter, start, help.length );
				start += help.length;
			}
		}
		
		public double[] getCurrentParameterValues( KindOfParameter kind ) throws Exception {
			if( parameter == null ) {
				fillParameters();
			}
			double[] res = parameter.clone();
			switch( kind ) {
				case PLUGIN:
					break;
				case LAST:
					break;
				case ZEROS:
					Arrays.fill( res, 0, score.length, 0 );
					break;
				default:
					throw new IllegalArgumentException( "Unknown kind of parameter" );
			}
			return res;
		}
		
		public void reset() throws Exception {
			prior.set(false, score);
		}
		
		protected KindOfParameter preoptimize( OptimizableFunction f ) throws Exception 
		{				
			/*TODO
			double workingSetSize = 0.1;
			DataSet[] original = f.getData();
			double[][] originalWeights = f.getSequenceWeights();
			NegativeDifferentiableFunction neg = new NegativeDifferentiableFunction(f);
			
			DataSet[] sampled = new DataSet[original.length];
			double[][] sampledWeights = new double[original.length][];		
			byte algo = (Byte) params.getParameterAt( 0 ).getValue();
			double t, linEps = (Double) params.getParameterAt( 2 ).getValue() * 1000d;
			StartDistanceForecaster sd = new ConstantStartDistance( (Double) params.getParameterAt( 3 ).getValue() );
			double[] parameter = f.getParameters( (KindOfParameter) params.getParameterAt( 5 ).getValue() );
//TODO eps = (Double) params.getParameterAt( 1 ).getValue() * 1000d
			TerminationCondition term = ((InstanceParameterSet<TerminationCondition>) params.getParameterAt( 1 ).getValue()).getInstance();//new SmallDifferenceOfFunctionEvaluationsCondition(eps);
			int it = 0;
			double[][] classified = ArrayHandler.clone(originalWeights);
			double[] help;
			DoubleList weights = new DoubleList();
			ArrayList<Sequence> seqList = new ArrayList<Sequence>();
			do {
				System.out.println("preoptimizing step " + it);
				
				//select new working set
				if( it == 0 ) {
					Pair<DataSet[], double[][]> p;
					for( int j = 0; j < sampled.length; j++ ) {
						p = original[j].partition( originalWeights[j], PartitionMethod.PARTITION_BY_WEIGHTS, 1d-workingSetSize, workingSetSize );
						sampled[j] = p.getFirstElement()[1];
						sampledWeights[j] = p.getSecondElement()[1];
					}
				} else {
					for( int clazz = 0; clazz < original.length; clazz++ ) {
						for( int seq = 0; seq < original[clazz].getNumberOfElements(); seq++ ) {
							classified[clazz][seq] = getLogProb( clazz, original[clazz].getElementAt(seq) );
						}
						//TODO
						help = classified[clazz].clone();
						Arrays.sort( help );
						t = help[(int)Math.ceil(help.length*workingSetSize)];
						System.out.println(clazz + ":\t" + t );
						seqList.clear();
						weights.clear();
						for( int seq = 0; seq < original[clazz].getNumberOfElements(); seq++ ) {
							if( classified[clazz][seq] <= t ) {
								seqList.add( original[clazz].getElementAt(seq) );
								weights.add( originalWeights[clazz][seq] );
							}
						}
						sampled[clazz] = new DataSet( "sampled", seqList.toArray( new Sequence[0] ) ); 
						sampledWeights[clazz] = weights.toArray();
					}
				}
				
				//set working set
				for( int clazz = 0; clazz < original.length; clazz++ ) {
					System.out.println( clazz + ":\t" + sampled[clazz].getNumberOfElements() );
				}
				f.setDataAndWeights(sampled, sampledWeights);
				
				//optimize
				Optimizer.optimize(algo, neg, parameter, term, linEps, sd, sostream );
				
				it++;
			} while( it < 1 );

			f.setDataAndWeights(original, originalWeights);
			System.out.println( "finished preoptimizing" );
			return KindOfParameter.LAST;
			/**/
//TODO
			return KindOfParameter.ZEROS;
		}
	}
}