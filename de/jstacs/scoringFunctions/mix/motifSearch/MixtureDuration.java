/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.scoringFunctions.mix.motifSearch;

import de.jstacs.NonParsableException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.data.DataSet;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;

/**
 * This class implements a mixture of {@link DurationScoringFunction}s.
 * 
 * @author Jens Keilwagen
 */
public class MixtureDuration extends DurationScoringFunction {

	private DurationScoringFunction[] function;
	private double[] hiddenParams, scores; 
	private double logNorm;
	private int[] paramRef, partDerOffset;
	private int starts;
	private IntList help;
	              
	
	private static double getESS( DurationScoringFunction... function ) {
		double ess = function[0].getESS(), e;
		boolean noESS = ess == 0;
		for( int i = 1; i < function.length; i++ ) {
			e = function[i].getESS();
			if( noESS ) {
				if( e > 0 ) {
					throw new IllegalArgumentException( "The ESS of duration " + i + " has to be zero." );
				}
			} else {
				ess += e;
			}
		}
		return ess;
	}
	
	/**
	 * The main constructor of a {@link MixtureDuration}.
	 * 
	 * @param starts the number of recommended starts
	 * @param function the {@link DurationScoringFunction}s for the components
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one element of <code>functions</code> could not
	 *             be cloned
	 * @throws IllegalArgumentException
	 *             if the <code>starts</code> is smaller than zero (0)
	 * @throws WrongAlphabetException
	 *             if at least one element of <code>function</code> has
	 *             an {@link de.jstacs.data.AlphabetContainer} that is not equal to those
	 *             of the other elements of <code>function</code>
	 */
	public MixtureDuration( int starts, DurationScoringFunction... function ) throws WrongAlphabetException, CloneNotSupportedException, IllegalArgumentException {
		super( function[0].getMin(), function[0].getMax(), getESS( function ) );
		if( starts > 0 ) {
			this.starts = starts;
		} else {
			throw new IllegalArgumentException( "The number of recommended starts should be positive." );
		}
		this.function = new DurationScoringFunction[function.length];
		scores = new double[function.length];
		hiddenParams = new double[function.length];
		paramRef = null;
		partDerOffset = new int[function.length];
		for( int i = 0; i < function.length; i++ ) {
			if( !alphabets.checkConsistency( function[i].getAlphabetContainer() ) ) {
				throw new WrongAlphabetException( "All durations have to have the same alphabet: Violated at position " + i );
			}
			this.function[i] = (DurationScoringFunction) function[i].clone();
		}
		logNorm = Normalisation.getLogSum( hiddenParams );
		setParamRef( false );
		help = new IntList();
	}
	
	private void setParamRef( boolean freeParams ) {
		if( paramRef == null || paramRef.length != function.length+2 ) {
			this.paramRef = new int[function.length+2];
		}
		int i = 0, n ;
		boolean unknown = false;
		for( ; i < function.length; i++ ) {
			n = function[i].getNumberOfParameters();
			unknown |= n < 0;
			paramRef[i+1] = paramRef[i] + n; 
		}
		if( unknown ) {
			paramRef[i+1] = -1;
		} else {
			paramRef[i+1] = paramRef[i] + scores.length - (freeParams?1:0);
		}
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link MixtureDuration} out of a {@link StringBuffer}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public MixtureDuration( StringBuffer source ) throws NonParsableException {
		super( source );
	}

	public MixtureDuration clone() throws CloneNotSupportedException {
		MixtureDuration clone = (MixtureDuration) super.clone();
		clone.function = ArrayHandler.clone( function );
		clone.scores = scores.clone();
		clone.paramRef = paramRef.clone();
		clone.partDerOffset = partDerOffset.clone();
		clone.hiddenParams = hiddenParams.clone();
		clone.help = new IntList();
		return clone;
	}
	
	@Override
	public void adjust( int[] length, double[] weight ) {
		double[][] assignedWeights = new double[function.length][];
		double[] stat = new double[hiddenParams.length];
		double all = 0;		
		for( int i = 0; i < function.length; i++ ) {
			function[i].adjust( length, weight ); //FIXME: problem if identical components
			hiddenParams[i] = 0;
			assignedWeights[i] = new double[weight.length];
			stat[i] = function[i].getESS();
			all += stat[i];
		}
		logNorm = Normalisation.getLogSum( hiddenParams );

				
		int[] values = new int[1];
		for( int l = 0; l < length.length; l++ ) {
			values[0] = length[l];
			for( int i = 0; i < function.length; i++ ) {
				scores[i] = hiddenParams[i] + function[i].getLogScore( values );
			}
			Normalisation.logSumNormalisation( scores );
			for( int i = 0; i < function.length; i++ ) {
				assignedWeights[i][l] = weight[l] * scores[i];
				stat[i] += weight[l] * scores[i];
			}
			all += weight[l];
		}

		for( int i = 0; i < function.length; i++ ) {
			function[i].adjust( length, assignedWeights[i] );
			hiddenParams[i] = Math.log( stat[i]/all );
		}
		logNorm = 0;
	}

	@Override
	public double getLogScore( int... values ) {
		for( int i = 0; i < function.length; i++ ) {
			scores[i] = hiddenParams[i] + function[i].getLogScore( values );
		}
		return Normalisation.getLogSum( scores ) - logNorm;
	}

	@Override
	public double getLogScoreAndPartialDerivation( IntList indices, DoubleList partialDer, int... values ) {
		int j, i = 0, o = partialDer.length();
		for( ; i < function.length; i++ ) {
			help.clear();
			scores[i] = hiddenParams[i] + function[i].getLogScoreAndPartialDerivation( help, partialDer, values );
			partDerOffset[i] = partialDer.length();	
			for( j = 0; j < help.length(); j++ ) {
				indices.add( paramRef[i] + help.get( j ) );
			}
		}
		double logScore = Normalisation.logSumNormalisation( scores );
		for( i = 0; i < function.length; i++ ) {
			partialDer.multiply( o, partDerOffset[i], scores[i] );
			o = partDerOffset[i];	
		}
		for( i = 0, j = paramRef[function.length]; j < paramRef[function.length+1]; j++, i++ ) {
			indices.add( j );
			partialDer.add( scores[i] - Math.exp( hiddenParams[i] - logNorm ) );	
		}
		return logScore - logNorm;
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		if( ess > 0 ) {
			for( int i = 0; i < function.length; i++ ) {
				function[i].addGradientOfLogPriorTerm( grad, paramRef[i] + start );
			}
			for( int j = 0, i = paramRef[function.length]; i < paramRef[function.length+1]; i++, j++ ) {
				grad[i+start] = function[j].getESS() - ess*Math.exp( hiddenParams[j] - logNorm );
			}
		}
	}

	public double getLogPriorTerm() {
		double lp = 0;
		if( ess > 0 ) {
			for( int i = 0; i < function.length; i++ ) {
				lp += function[i].getLogPriorTerm()
					+ function[i].getESS() * hiddenParams[i];
			}
			lp -= ess * logNorm;
		}
		return lp;
	}

	public double[] getCurrentParameterValues() throws Exception {
		int n = getNumberOfParameters(); 
		if( n > 0 ) {
			double[] params = new double[n], current;
			for( int i = 0; i < function.length; i++ ) {
				current = function[i].getCurrentParameterValues();
				System.arraycopy( current, 0, params, paramRef[i], current.length );
			}
			for( int j = 0, i = paramRef[function.length]; i < paramRef[function.length+1]; i++, j++ ) {
				params[i] = hiddenParams[j];
			}
			return params;
		} else {
			throw new RuntimeException();
		}
	}
	
	public String getInstanceName() {
		String name = "mixture duration(" + function[0].getInstanceName();
		for( int i = 1; i < function.length; i++ ) {
			name += ", " + function[i].getInstanceName();
		}
		return name + ")";
	}

	public int getNumberOfParameters() {
		return paramRef[function.length+1];
	}

	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		int i = 0, j;
		double w = 1, all = 0;
		for( ; i < function.length; i++ ) {
			hiddenParams[i] = function[i].getESS();
			all += hiddenParams[i];
		}
		DirichletMRGParams params = new DirichletMRGParams( hiddenParams );
		double[] help, current = new double[function.length];
		if( weights == null ) {
			weights = new double[data.length][];
		}
		help = weights[index];
		double[][] componentWeights = new double[function.length][data[index].getNumberOfElements()];
		for( j = 0; j < componentWeights[0].length; j++ ) {
			DirichletMRG.DEFAULT_INSTANCE.generate( current, 0, function.length, params );
			if( help != null ) {
				w = help[j];
			}
			all += w;
			//System.out.println( w + "\t"  + Arrays.toString( current ) );
			for( i = 0; i < function.length; i++ ) {
				componentWeights[i][j] = w*current[i];
				hiddenParams[i] += componentWeights[i][j];
			}
		}
		for( i = 0; i < function.length; i++ ) {
			weights[index] = componentWeights[i];
			function[i].initializeFunction( index, freeParams, data, weights );
			//System.out.println( hiddenParams[i] / all );
			hiddenParams[i] = Math.log( hiddenParams[i] / all );
		}
		logNorm = 0;
		weights[index] = help;
		setParamRef( freeParams );
	}
	
	@Override
	public void initializeUniformly() {
		for( int i = 0; i < function.length; i++ ) {
			function[i].initializeUniformly();
			hiddenParams[i] = 0;
		}
		logNorm = Normalisation.getLogSum( hiddenParams );
		setParamRef( false ); 
	}

	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		boolean noPrior = ess == 0;
		double[] hyper = scores.clone();
		for( int i = 0; i < function.length; i++ ) {
			function[i].initializeFunctionRandomly( freeParams );
			hyper[i] = noPrior ? 1 : function[i].getESS();
		}
		DirichletMRG.DEFAULT_INSTANCE.generate( hiddenParams, 0, function.length, new DirichletMRGParams( hyper ) );
		for( int i = 0; i < function.length; i++ ) {
			hiddenParams[i] = Math.log( hiddenParams[i] );
		}
		logNorm = 0;
		setParamRef( freeParams );
	}

	public boolean isInitialized() {
		int i = 0;
		while( i < function.length && function[i].isInitialized() ) {
			i++;
		}
		return i == function.length;
	}
	
	public boolean isNormalized() {
		return true;
	}

	public void setParameters( double[] params, int start ) {
		for( int i = 0; i < function.length; i++ ) {
			function[i].setParameters( params, start + paramRef[i] );
		}
		for( int j = 0, i = paramRef[function.length]; i < paramRef[function.length+1]; i++, j++ ) {
			hiddenParams[j] = params[start+i];
		}
		logNorm = Normalisation.getLogSum( hiddenParams );
	}

	@Override
	protected String getRNotation( String distributionName ) {
		String r = "", sum = null;
		for( int i = 0; i < function.length; i++ ) {
			r += function[i].getRNotation( distributionName + i ) + "\n";
			if( sum == null ) {
				sum = distributionName + " = " + Math.exp( hiddenParams[i] - logNorm ) + " * " + distributionName + i;  
			} else {
				sum += " + " + Math.exp( hiddenParams[i] - logNorm ) + " * " + distributionName + i;
			}
		}
		return r + sum  + ";";
	}
	
	public void modify( int delta ) {
		if( delta != 0 ) {
			super.modify( delta );
			for( int i =0; i < function.length; i++ ) {
				function[i].modify( delta );
			}
		}
	}
	
	public int getNumberOfRecommendedStarts() {
		return starts;
	}
		
	protected void fromXML( StringBuffer rep ) throws NonParsableException
	{
		StringBuffer xml = XMLParser.extractForTag( rep, XML_TAG );
		super.fromXML(xml);
		function = XMLParser.extractObjectForTags( xml, "components", DurationScoringFunction[].class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		hiddenParams = XMLParser.extractObjectForTags( xml, "hiddenParams", double[].class );
		starts = XMLParser.extractObjectForTags( xml, "starts", int.class );

		scores = new double[function.length];
		paramRef = null;
		partDerOffset = new int[function.length];
		logNorm = Normalisation.getLogSum( hiddenParams );
		setParamRef( false );
		help = new IntList();
	}
	
	public StringBuffer toXML()
	{
		StringBuffer xml = super.toXML();
		XMLParser.appendObjectWithTags( xml, function, "components" );
		XMLParser.appendObjectWithTags( xml, hiddenParams, "hiddenParams" );
		XMLParser.appendObjectWithTags( xml, starts, "starts" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}
	
	private static String XML_TAG = "MixtureDuration";
}
