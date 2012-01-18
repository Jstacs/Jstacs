package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions;

import java.text.NumberFormat;
import java.util.Arrays;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DiMRGParams;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jstacs.utils.random.FastDirichletMRGParams;

/**
 * This class implements a mixture of {@link Emission}s.
 * 
 * @author Jens Keilwagen
 */
public final class MixtureEmission implements Emission {

	private Emission[] emission;
	private double[] hyper;
	private double[] logProb;
	
	private double[] help, statistic;
	
	/**
	 * The main constructor creating a {@link MixtureEmission} from a set of emissions.
	 * 
	 * @param emission the individual emissions
	 * @param hyperParameters the hyper parameters for each component
	 * 
	 * @throws CloneNotSupportedException if the emission could not be cloned
	 */
	public MixtureEmission( Emission[] emission, double[] hyperParameters ) throws CloneNotSupportedException {
		if( hyperParameters != null && emission.length != hyperParameters.length ) {
			throw new IllegalArgumentException( "The number of emissions and the number of hyper-parameters has to be equal." );
		}
		AlphabetContainer con = emission[0].getAlphabetContainer();
		for( int e = 1; e < emission.length; e++ ) {
			if( !con.checkConsistency( emission[e].getAlphabetContainer() ) ) {
				throw new IllegalArgumentException( "The emissions must work on the same AlphabetContainer." );
			}
			if( hyperParameters != null && hyperParameters[e] < 0 ) {
				throw new IllegalArgumentException( "The hyper-parameters have to be non-negative." );
			}
		}
		this.emission = ArrayHandler.clone( emission );
		hyper = hyperParameters== null ? new double[emission.length] : hyperParameters.clone();
		logProb = new double[emission.length];
		Arrays.fill( logProb, -Math.log( emission.length ) );
		init();
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link MixtureEmission} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MixtureEmission} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public MixtureEmission( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, XML_TAG );
		emission = (Emission[]) XMLParser.extractObjectForTags( xml, "emissions" );
		hyper = (double[]) XMLParser.extractObjectForTags( xml, "hyper" );
		logProb = (double[]) XMLParser.extractObjectForTags( xml, "logProb" );
		init();
	}
	
	private static final String XML_TAG = "MixtureEmission";
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, emission, "emissions" );
		XMLParser.appendObjectWithTags( xml, hyper, "hyper" );
		XMLParser.appendObjectWithTags( xml, logProb, "logProb" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}
	
	private void init() {
		help = new double[emission.length];
		statistic = new double[emission.length];
	}
	
	public MixtureEmission clone() throws CloneNotSupportedException {
		MixtureEmission clone = (MixtureEmission) super.clone();
		clone.emission = ArrayHandler.clone( emission );
		clone.help = clone.help;
		clone.hyper = hyper.clone();
		clone.logProb = logProb.clone();
		clone.statistic = statistic.clone();
		return clone;
	}

	public void joinStatistics(Emission... emissions){
		Emission[] temp = new Emission[emissions.length];
		for(int i=0;i<this.emission.length;i++){
			for(int j=0;j<emissions.length;j++){
				temp[j] = ((MixtureEmission)emissions[j]).emission[j];
			}
			this.emission[i].joinStatistics( temp );
		}
		for(int i=0;i<emissions.length;i++){
			if(emissions[i] != this){
				for(int j=0;j<statistic.length;j++){
					statistic[j] += ((MixtureEmission)emissions[i]).statistic[j];
				}
			}
		}
		for(int i=0;i<emissions.length;i++){
			if(emissions[i] != this){
				System.arraycopy( this.statistic, 0, ((MixtureEmission)emissions[i]).statistic, 0, this.statistic.length );
			}
		}
	}
	
	@Override
	public void addToStatistic( boolean forward, int startPos, int endPos, double weight, Sequence seq ) throws OperationNotSupportedException {
		for( int e = 0; e < emission.length; e++ ) {
			help[e] = logProb[e] - emission[e].getLogProbFor( forward, startPos, endPos, seq );
		}
		Normalisation.logSumNormalisation( help );
		for( int e = 0; e < emission.length; e++ ) {
			emission[e].addToStatistic( forward, startPos, endPos, weight*help[e], seq );
			statistic[e] += help[e];
		}
	}

	@Override
	public void estimateFromStatistic() {
		double sum = 0;
		for( int e = 0; e < emission.length; e++ ) {
			emission[e].estimateFromStatistic();
			statistic[e] += hyper[e];
			sum += statistic[e];
		}
		sum = Math.log( sum );
		for( int e = 0; e < emission.length; e++ ) {
			logProb[e] = Math.log( statistic[e] ) - sum;
		}
	}

	@Override
	public AlphabetContainer getAlphabetContainer() {
		return emission[0].getAlphabetContainer();
	}

	@Override
	public double getLogPriorTerm() {
		double lp = 0;
		for( int e = 0; e < emission.length; e++ ) {
			lp += emission[e].getLogPriorTerm();
			lp += logProb[e] * hyper[e]; 
		}
		return lp;
	}

	@Override
	public double getLogProbFor( boolean forward, int startPos, int endPos, Sequence seq ) throws OperationNotSupportedException {
		for( int e = 0; e < emission.length; e++ ) {
			help[e] = logProb[e] - emission[e].getLogProbFor( forward, startPos, endPos, seq );
		}
		return Normalisation.getLogSum( 0, emission.length, help );
	}

	@Override
	public void initializeFunctionRandomly() {
		int zero = 0;
		for( int e = 0; e < emission.length; e++ ) {
			if( hyper[e] == 0 ) {
				zero++;
			}
		}
		DiMRGParams p;
		if( zero == emission.length ) {
			p = new FastDirichletMRGParams( 1 );
		} else if ( zero == 0 ) {
			p = new DirichletMRGParams( hyper );
		} else {
			throw new IllegalArgumentException();
		}
		
		DirichletMRG.DEFAULT_INSTANCE.generateLog( logProb, 0, emission.length, p );
		for( int e = 0; e < emission.length; e++ ) {
			emission[e].initializeFunctionRandomly();
		}
	}

	@Override
	public void resetStatistic() {
		for( int e = 0; e < emission.length; e++ ) {
			emission[e].resetStatistic();
		}
	}

	@Override
	public String getNodeShape(boolean forward) {
		String res = "";
		if( getAlphabetContainer().isReverseComplementable() ) {
			res += "\"house\", orientation=";
			if ( forward ){
				res+= "-";
			}
			res +="90";
		} else {
			res += "\"box\"";
		}
		return res;
	}

	@Override
	public String getNodeLabel( double weight, String name, NumberFormat nf ) {
		return "\""+name+"\"";
	}
	
	
	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		if( !t.getClass().equals( getClass() ) ) {
			throw new IllegalArgumentException( "The emissions are not comparable." );
		}
		MixtureEmission tt = (MixtureEmission) t;
		for( int i = 0; i < emission.length; i++ ) {
			emission[i].setParameters( tt.emission[i] );
		}		
	}
}
