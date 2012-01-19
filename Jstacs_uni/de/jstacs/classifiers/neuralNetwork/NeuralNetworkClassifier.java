package de.jstacs.classifiers.neuralNetwork;

import java.util.Random;

import de.jstacs.NotTrainedException;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.classifiers.neuralNetwork.activationFunctions.ActivationFunction;
import de.jstacs.classifiers.neuralNetwork.neurons.InnerNeuron;
import de.jstacs.classifiers.neuralNetwork.neurons.InputNeuron;
import de.jstacs.classifiers.neuralNetwork.neurons.MSEOutputNeuron;
import de.jstacs.classifiers.neuralNetwork.neurons.Neuron;
import de.jstacs.classifiers.neuralNetwork.stepSizeAdaption.StepSizeAdaption;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.Time;


public class NeuralNetworkClassifier extends AbstractScoreBasedClassifier {

	private static Random r = new Random();
	
	public enum Learning{
		BATCH,
		STOCHASTIC
	};
	
	private Neuron[][] neurons;
	private boolean isInitialized;
	private int[][] freeIndexes;
	private int lastIndex;
	private Learning learning;
	private StepSizeAdaption stepSize;
	private ActivationFunction activationFunction;
	private double initStepSize;
	private TerminationCondition termCond;
	
	public NeuralNetworkClassifier( AlphabetContainer abc, int[] numNeurons, ActivationFunction activationFunction, Learning learning, StepSizeAdaption stepSize, double initStepSize, TerminationCondition termCond ) {
		super( abc, numNeurons[numNeurons.length-1] );
		neurons = new Neuron[numNeurons.length][];
		for(int i=0;i<numNeurons.length;i++){
			neurons[i] = new Neuron[numNeurons[i]];
			for(int j=0;j<numNeurons[i];j++){
				if(i==0){
					neurons[i][j] = new InputNeuron( j );
				}else if(i==numNeurons.length-1){
					neurons[i][j] = new MSEOutputNeuron( activationFunction, j, neurons[i-1] );
				}else{
					neurons[i][j] = new InnerNeuron( activationFunction, j, neurons[i-1] );
				}
			}
		}
		this.activationFunction = activationFunction;
		this.learning = learning;
		this.stepSize = stepSize;
		this.initStepSize = initStepSize;
		this.termCond = termCond;
		isInitialized = false;
	}

	public NeuralNetworkClassifier( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	protected double getScore( Sequence seq, int i, boolean check ) throws IllegalArgumentException, NotTrainedException, Exception {
		return neurons[neurons.length-1][i].getOutput( seq );
	}

	@Override
	public String getInstanceName() {
		return "Neural network classifier";
	}

	@Override
	public CategoricalResult[] getClassifierAnnotation() {
		return null;
	}

	@Override
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		return null;
	}

	@Override
	public boolean isInitialized() {
		return isInitialized;
	}
	
	protected int[][] getIndexes( ){
		if(learning == Learning.BATCH){
			return freeIndexes;
		}else{
			int d = r.nextInt( lastIndex );
			int[][] temp = new int[][]{freeIndexes[d]};
			lastIndex--;
			freeIndexes[d] = freeIndexes[lastIndex];
			freeIndexes[lastIndex] = temp[0];
			if(lastIndex == 0){
				lastIndex = freeIndexes.length;
			}
			return temp;
		}
	}
	
	protected void prepareIndexes( DataSet[] data ){
		int num = 0;
		for(int i=0;i<data.length;i++){
			num += data[i].getNumberOfElements();
		}
		freeIndexes = new int[num][2];
		for(int i=0,j=0;i<data.length;i++){
			for(int k=0;k<data[i].getNumberOfElements();k++,j++){
				freeIndexes[j][0] = i;
				freeIndexes[j][1] = k;
			}
		}
		lastIndex = freeIndexes.length;
	}

	@Override
	public void train( DataSet[] data, double[][] weights ) throws Exception {
		
		for(int i=0;i<neurons.length;i++){
			for(int j=0;j<neurons[i].length;j++){
				neurons[i][j].initializeRandomly();
			}
		}
		isInitialized = true;
		
		double[][] desiredOutputs = new double[data.length][data.length];
		for(int i=0;i<desiredOutputs.length;i++){
			for(int j=0;j<desiredOutputs[i].length;j++){
				desiredOutputs[i][j] = i==j ? activationFunction.getPositiveValue() : activationFunction.getNegativeValue();
			}
		}
		
		prepareIndexes( data );
		
		int epoch = 0;
		int iteration = 1;
		int num = 0;
		double totalError = Double.POSITIVE_INFINITY;
		double lastTotalError = totalError;
		double currStepSize = initStepSize;

		Time time = Time.getTimeInstance( null );
		
		do{
			if(epoch % 1000 == 0){
				System.out.println(epoch+"\t"+totalError+"\t"+time.getElapsedTime());
			}
			num = freeIndexes.length;
			epoch++;
			lastTotalError = totalError;
			totalError = 0;

			while(num>0){

				int[][] idxs = getIndexes( );

				for(int n=0;n<idxs.length;n++){

					for(int i=0;i<neurons.length;i++){
						for(int j=0;j<neurons[i].length;j++){
							neurons[i][j].reset();
						}
					}

					Sequence input = data[idxs[n][0]].getElementAt( idxs[n][1] );
					double weight = weights == null || weights[idxs[n][0]] == null ? 1 : weights[idxs[n][0]][idxs[n][1]];

					for(int i=0;i<neurons[neurons.length-1].length;i++){
						double out = neurons[neurons.length-1][i].getOutput( input );
						/*if(Math.abs( desiredOutputs[idxs[n][0]][i]-out ) > 1E-6){
						System.out.println(n+": "+desiredOutputs[idxs[n][0]][i]+" "+out+" "+((InnerNeuron)neurons[neurons.length-1][i]).h);
					}*/
						totalError += 0.5*weight*(desiredOutputs[idxs[n][0]][i]-out)*(desiredOutputs[idxs[n][0]][i]-out);

					}
					for(int i=0;i<neurons[0].length;i++){
						neurons[0][i].getError( input, weight, desiredOutputs[idxs[n][0]] );
						//System.out.println(n+": "+((InnerNeuron)neurons[neurons.length-1][i]).error);
					}

				}

				for(int i=0;i<neurons.length;i++){
					for(int j=0;j<neurons[i].length;j++){
						neurons[i][j].adaptWeights( currStepSize );
					}
				}

				iteration++;
				num -= idxs.length;


				currStepSize = stepSize.getStepSize( initStepSize, currStepSize, iteration, epoch );

			}
			
			
		}while(termCond.doNextIteration( epoch, lastTotalError, totalError, null, null, currStepSize, time ));
		
		System.out.println(epoch+"\t"+totalError+"\t"+time.getElapsedTime());
		
		freeIndexes = null;
		
	}

	@Override
	protected String getXMLTag() {
		return getClass().getSimpleName();
	}

	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, activationFunction, "activationFunction" );
		XMLParser.appendObjectWithTags( xml, initStepSize, "initStepSize" );
		XMLParser.appendObjectWithTags( xml, isInitialized, "isInitialized" );
		XMLParser.appendObjectWithTags( xml, learning, "learning" );
		XMLParser.appendObjectWithTags( xml, neurons, "neurons" );
		XMLParser.appendObjectWithTags( xml, stepSize, "stepSize" );
		XMLParser.appendObjectWithTags( xml, termCond, "termCond" );
		return xml;
	}

	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		activationFunction = XMLParser.extractObjectForTags( xml, "activationFunction", ActivationFunction.class );
		initStepSize = XMLParser.extractObjectForTags( xml, "initStepSize", double.class );
		isInitialized = XMLParser.extractObjectForTags( xml, "isInitialized", boolean.class );
		learning = XMLParser.extractObjectForTags( xml, "learning", Learning.class );
		neurons = XMLParser.extractObjectForTags( xml, "neurons", Neuron[][].class );
		stepSize = XMLParser.extractObjectForTags( xml, "stepSize", StepSizeAdaption.class );
		termCond = XMLParser.extractObjectForTags( xml, "termCond", TerminationCondition.class );
		
		for(int i=1;i<neurons.length;i++){
			for(int j=0;j<neurons[i].length;j++){
				((InnerNeuron)neurons[i][j]).setPredecessors( neurons[i-1] );
			}
		}
		
	}
	
	
	

}
