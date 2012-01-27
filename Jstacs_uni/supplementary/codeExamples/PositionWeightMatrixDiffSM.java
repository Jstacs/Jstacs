package supplementary.codeExamples;

import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;


public class PositionWeightMatrixDiffSM extends AbstractDifferentiableStatisticalModel {

	private double[][] parameters;// array for the parameters of the PWM in natural parameterization
	private double ess;
	private boolean isInitialized;
	private Double norm;
	
	public PositionWeightMatrixDiffSM( AlphabetContainer alphabets, int length, double ess ) throws IllegalArgumentException {
		super( alphabets, length );
		if(!alphabets.isSimple() || !alphabets.isDiscrete()){
			throw new IllegalArgumentException( "This PWM can handle only discrete alphabets with the same alphabet at each position." );
		}
		//create parameter-array
		this.parameters = new double[length][(int)alphabets.getAlphabetLengthAt( 0 )];
		this.ess = ess;
		this.isInitialized = false;
	}

	/**
	 * @param xml
	 * @throws NonParsableException
	 */
	public PositionWeightMatrixDiffSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		return parameters[0].length;
	}

	@Override
	public double getLogNormalizationConstant() {
		if(this.norm == null){
			norm = 0.0;
			for(int i=0;i<parameters.length;i++){
				norm += Normalisation.getLogSum( parameters[i] );
			}
		}
		return norm;
	}

	@Override
	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception {
		if(norm == null){
			getLogNormalizationConstant();
		}
		int row = parameterIndex;
		int column = 0;
		while(row > parameters[column].length){
			row -= parameters[column].length;
			column++;
		}
		return norm - Normalisation.getLogSum( parameters[column] ) + parameters[column][row];
	}

	@Override
	public double getLogPriorTerm() {
		double logPrior = 0;
		for(int i=0;i<parameters.length;i++){
			for(int j=0;j<parameters[i].length;j++){
				logPrior += ess/alphabets.getAlphabetLengthAt( 0 ) * parameters[i][j];
			}
		}
		return logPrior;
	}

	@Override
	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		for(int i=0;i<parameters.length;i++){
			for(int j=0;j<parameters[i].length;j++,start++){
				grad[start] = ess/alphabets.getAlphabetLengthAt( 0 );
			}
		}
	}

	@Override
	public double getESS() {
		return ess;
	}

	@Override
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		if(!data[index].getAlphabetContainer().checkConsistency( alphabets ) || 
				data[index].getElementLength() != length){
			throw new IllegalArgumentException( "Alphabet or length to not match." );
		}
		for(int i=0;i<parameters.length;i++){
			Arrays.fill( parameters[i], ess/alphabets.getAlphabetLengthAt( 0 ) );
		}
		for(int i=0;i<data[index].getNumberOfElements();i++){
			Sequence seq = data[index].getElementAt( i );
			for(int j=0;j<seq.getLength();j++){
				parameters[j][ seq.discreteVal( j ) ] += weights[index][i];
			}
		}
		for(int i=0;i<parameters.length;i++){
			Normalisation.sumNormalisation( parameters[i] );
			for(int j=0;j<parameters[i].length;j++){
				parameters[i][j] = Math.log( parameters[i][j] );
			}
		}
		norm = null;
		isInitialized = true;
	}

	@Override
	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		int al = (int)alphabets.getAlphabetLengthAt( 0 );
		DirichletMRGParams pars = new DirichletMRGParams( ess/al, al );
		for(int i=0;i<parameters.length;i++){
			parameters[i] = DirichletMRG.DEFAULT_INSTANCE.generate( al, pars );
		}
		norm = null;
		isInitialized = true;
	}
	

	@Override
	public double getLogScoreFor( Sequence seq, int start ) {
		double score = 0.0;
		for(int i=0;i<parameters.length;i++){
			score += parameters[i][ seq.discreteVal( i+start ) ];
		}
		return score;
	}

	@Override
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		double score = 0.0;
		int off = 0;
		for(int i=0;i<parameters.length;i++){
			int v = seq.discreteVal( i+start );
			score += parameters[i][ v ];
			indices.add( off + v );
			partialDer.add( 1 );
			off += parameters[i].length;
		}
		return score;
	}

	@Override
	public int getNumberOfParameters() {
		int num = 0;
		for(int i=0;i<parameters.length;i++){
			num += parameters[i].length;
		}
		return num;
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] pars = new double[getNumberOfParameters()];
		for(int i=0,k=0;i<parameters.length;i++){
			for(int j=0;j<parameters[i].length;j++,k++){
				pars[k] = parameters[i][j];
			}
		}
		return pars;
	}

	@Override
	public void setParameters( double[] params, int start ) {
		for(int i=0;i<parameters.length;i++){
			for(int j=0;j<parameters[i].length;j++,start++){
				parameters[i][j] = params[start];
			}
		}
		norm = null;
	}

	@Override
	public String getInstanceName() {
		return "Position weight matrix";
	}


	@Override
	public boolean isInitialized() {
		return isInitialized;
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, parameters, "parameters" );
		XMLParser.appendObjectWithTags( xml, isInitialized, "isInitialized" );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
		XMLParser.addTags( xml, "PWM" );
		return xml;
	}

	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, "PWM" );
		alphabets = (AlphabetContainer)XMLParser.extractObjectForTags( xml, "alphabets" );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		parameters = (double[][])XMLParser.extractObjectForTags( xml, "parameters" );
		isInitialized = XMLParser.extractObjectForTags( xml, "isInitialized", boolean.class );
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
	}

}
