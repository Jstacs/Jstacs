/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.scoringFunctions;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.scoringFunctions.homogeneous.HomogeneousScoringFunction;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class enables the user to model parts of a sequence independent of each
 * other. For instance, the first part of the sequence is modeled by the first
 * {@link NormalizableScoringFunction} and has the length of the first
 * {@link NormalizableScoringFunction}, the second part starts directly after
 * the first part, is modeled by the second {@link NormalizableScoringFunction}
 * ... etc. It is also possible to use a {@link NormalizableScoringFunction} for
 * more than one sequence part and in both orientations (if possible).
 * 
 * <br><br>
 * 
 * It is important to set the equivalent sample size (ESS) of each instance carefully, i.e., corresponding to the ESS of the parts.  
 * 
 * @author Jens Keilwagen
 */
public class BasicIndependentProductScoringFunction extends AbstractScoringFunction {
	
	protected ScoringFunction[] score;

	protected int[] index, start, partialLength, params;
	
	protected boolean[] reverse, isVariable;
	
	private IntList partIList;
	
	private boolean plugIn;

	private final static AlphabetContainer getAlphabetContainer( ScoringFunction[] functions, int[] index, int length[], boolean[] reverse ) throws IllegalArgumentException, WrongAlphabetException {
		boolean direct = index == null;
		int len = direct ? functions.length : index.length;
		AlphabetContainer[] cons = new AlphabetContainer[len];
		int[] lengths = new int[len];
		for( int j, i = 0; i < len; i++ ) {
			if( direct ) {
				j = i;
			} else {
				j = index[i];
			}
			cons[i] = functions[j].getAlphabetContainer();
			lengths[i] = length[i];
			if( reverse != null && reverse[i] == true && !cons[i].isReverseComplementable() ) {
				throw new WrongAlphabetException( "The AlpabetContainer is not reverse complementable." );
			}
		}
		return new AlphabetContainer( cons, lengths );
	}

	private final static int sum( int[] length ) throws IllegalArgumentException {
		int res = 0, i = 0;
		while( i < length.length && length[i] > 0 ) {
			res += length[i++];
		}
		if( i != length.length ) {
			throw new IllegalArgumentException( "The length with index " + i + " is 0." );
		}
		return res;
	}

	protected final static int[] getLengthArray( ScoringFunction... functions ) throws IllegalArgumentException {
		int i = 0;
		int[] res = new int[functions.length];
		while( i < functions.length && functions[i].getLength() > 0 ) {
			res[i] = functions[i].getLength();
			i++;
		}
		if( i != functions.length ) {
			throw new IllegalArgumentException( "The ScoringFunction with index " + i + " has a length 0." );
		}
		return res;
	}

	/**
	 * This constructor creates an instance of an
	 * {@link BasicIndependentProductScoringFunction} from a given series of
	 * independent {@link ScoringFunction}s. The length that is
	 * modeled by each component is determined by
	 * {@link ScoringFunction#getLength()}. So the length should not be 0.
	 * 
	 * @param plugIn whether to use plugIn parameters for the parts, otherwise the last parameters are used for parts that are instance of {@link HomogeneousScoringFunction}
	 * @param functions
	 *            the components, i.e. the given series of independent
	 *            {@link NormalizableScoringFunction}s
	 *             
	 * @throws CloneNotSupportedException
	 *             if at least one element of <code>functions</code> could not
	 *             be cloned
	 * @throws WrongAlphabetException
	 *             if the user tries to use an alphabet for a reverse complement that can not be used for a reverse complement.
	 *             
	 * @see BasicIndependentProductScoringFunction#BasicIndependentProductScoringFunction(boolean, ScoringFunction[], int[])
	 */
	public BasicIndependentProductScoringFunction( double ess, boolean plugIn, ScoringFunction... functions ) throws CloneNotSupportedException, WrongAlphabetException {
		this( plugIn, functions, getLengthArray( functions ) );
	}

	/**
	 * This constructor creates an instance of an
	 * {@link BasicIndependentProductScoringFunction} from given series of
	 * independent {@link ScoringFunction}s and lengths.
	 * 
	 * @param plugIn whether to use plugIn parameters for the parts, otherwise the last parameters are used for parts that are instance of {@link HomogeneousScoringFunction}
	 * @param functions
	 *            the components, i.e. the given series of independent
	 *            {@link ScoringFunction}s
	 * @param length
	 *            the lengths, one for each component
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one component could not be cloned
	 * @throws WrongAlphabetException
	 *             if the user tries to use an alphabet for a reverse complement that can not be used for a reverse complement.
	 * 
	 * @see BasicIndependentProductScoringFunction#BasicIndependentProductScoringFunction(boolean, ScoringFunction[], int[], int[], boolean[])
	 */
	public BasicIndependentProductScoringFunction( boolean plugIn, ScoringFunction[] functions, int[] length ) throws CloneNotSupportedException, WrongAlphabetException {
		this( plugIn, functions, null, length, null );
	}

	/**
	 * This is the main constructor.
	 *
	 * @param plugIn whether to use plugIn parameters for the parts, otherwise the last parameters are used for parts that are instance of {@link HomogeneousScoringFunction}
	 * @param functions the {@link ScoringFunction}
	 * @param index the index of the {@link ScoringFunction} at each part
	 * @param length the length of each part
	 * @param reverse a switch whether to use it directly or the reverse complementary strand
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one component could not be cloned
	 * @throws WrongAlphabetException
	 *             if the user tries to use an alphabet for a reverse complement that can not be used for a reverse complement. 
	 * 
	 */
	public BasicIndependentProductScoringFunction( boolean plugIn, ScoringFunction[] functions, int[] index, int[] length, boolean[] reverse ) throws CloneNotSupportedException, WrongAlphabetException {
		super( getAlphabetContainer( functions, index, length, reverse ), sum( length ) );
		this.plugIn = plugIn;
		score = ArrayHandler.clone( functions );
		set( index, length, reverse );
		setParamsStarts();
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link BasicIndependentProductScoringFunction} out of a
	 * {@link StringBuffer} as returned by {@link #toXML()}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public BasicIndependentProductScoringFunction( StringBuffer source ) throws NonParsableException {
		super( source );
	}

	private void set( int[] index, int[] length, boolean[] reverse ) throws IllegalArgumentException {
		int[] used = new int[score.length];
		boolean direct = index == null;
		int oldStart = 0, len = direct ? score.length : index.length;
		start = new int[len];
		partialLength = new int[len];
		isVariable = new boolean[score.length];
		this.index = new int[len];
		this.reverse = new boolean[len];
		for( int i = 0; i < len; i++ ) {
			if( direct ) {
				this.index[i] = i;
			} else {
				if( 0 <= index[i] && index[i] < score.length ) {
					this.index[i] = index[i]; 
				} else {
					throw new IndexOutOfBoundsException( "index " + index[i] );
				}
			}
			used[this.index[i]]++;
			if( reverse == null ) {
				this.reverse[i] = false;
			} else {
				this.reverse[i] = reverse[i];
			}
			
			start[i] = oldStart;
			partialLength[i] = length[i];
			isVariable[this.index[i]] = score[this.index[i]] instanceof VariableLengthScoringFunction;
			if( !isVariable[this.index[i]] && score[this.index[i]].getLength() != partialLength[i] ) {
				throw new IllegalArgumentException( "Could not use length " + partialLength[i] + " at part " + i + " for ScoringFunction with index " + this.index[i] + "." );
			}
			oldStart += partialLength[i];
		}
		for( int i = 0; i < used.length; i++ ) {
			if( used[i] == 0 ) {
				throw new IllegalArgumentException( "The ScoringFunction with index " + i + " is never used." );
			}
		}
		partIList = new IntList();
	}

	protected void setParamsStarts() {
		if( params == null ) {
			params = new int[score.length + 1];
		}
		for( int n, i = 0; i < score.length; i++ ) {
			n = score[i].getNumberOfParameters();
			if( n == UNKNOWN ) {
				params = null;
				break;
			} else {
				params[i + 1] = params[i] + n;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#clone()
	 */
	public BasicIndependentProductScoringFunction clone() throws CloneNotSupportedException {
		BasicIndependentProductScoringFunction clone = (BasicIndependentProductScoringFunction)super.clone();
		clone.score = ArrayHandler.clone( score );
		clone.set( index, partialLength, reverse );
		clone.params = null;
		clone.setParamsStarts();
		return clone;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int, boolean, de.jstacs.data.Sample[], double[][])
	 */
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		DataSet[] part = new DataSet[data.length];
		double[][] help;
		for( int a, i = 0; i < score.length; i++ ) {
			if( plugIn || !(score[i] instanceof HomogeneousScoringFunction) ) {
				a = extractSequenceParts( i, data, part );
				help = a==1 ? weights : extractWeights( a, weights );
				score[i].initializeFunction( index, freeParams, part, help );
			}
		}		
		setParamsStarts();
	}
	
	/**
	 * This method extracts the corresponding {@link Sequence} parts for a specific {@link ScoringFunction}.
	 * 
	 * @param scoringFunctionIndex the index of the {@link ScoringFunction}
	 * @param data the original data
	 * @param result an array for the resulting {@link DataSet}s of {@link Sequence}s; has to have same length as <code>data</code>
	 * 
	 * @return the number how often the {@link ScoringFunction} was used
	 * 
	 * @throws Exception if the Sample can not be created
	 */
	public int extractSequenceParts( int scoringFunctionIndex, DataSet[] data, DataSet[] result ) throws Exception {
		DataSet current;
		Arrays.fill( result, null );
		int used = 0;
		for( int n, j, k = 0; k < index.length; k++ ) {
			if( index[k] == scoringFunctionIndex ) {
				used++;
				for( j = 0; j < data.length; j++ ) {
					if( data[j] != null ) {
						current = data[j].getInfixDataSet( start[k], partialLength[k] );
						if( reverse[k] ) {
							Sequence[] seq = new Sequence[current.getNumberOfElements()];
							for( n = 0; n < seq.length; n++ ) {
								seq[n] = current.getElementAt( n ).reverseComplement();
							}
							current = new DataSet( "reverse complement of \"" + current.getAnnotation() +"\"", seq );
						}
						if( result[j] == null ) {
							result[j] = current;
						} else {
							result[j] = DataSet.union( result[j], current );
						}
					}
				}
			}
		}
		return used;
	}
	
	/**
	 * This method creates the weights for {@link BasicIndependentProductScoringFunction#extractSequenceParts(int, DataSet[], DataSet[])}.
	 * 
	 * @param number the number how often the weights should be copied after each other.
	 * @param weights the original weights
	 * 
	 * @return the new weights (might be <code>null</code>)
	 * 
	 * @see #extractSequenceParts(int, DataSet[], DataSet[])
	 */
	public double[][] extractWeights( int number, double[][] weights ) {
		double[][] res;
		if( number == 1 || weights == null ) {
			res = weights;
		} else {
			res = new double[weights.length][];
			for( int n, j, k = 0; k < res.length; k++ ) {
				n = weights[k].length;
				res[k] = new double[number*n];
				for( j = 0; j < number; j++ ) {
					System.arraycopy( weights[k], 0, res[k], j*n, n );
				}
			}
		}
		return res;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#fromXML(java.lang.StringBuffer)
	 */
	protected void fromXML( StringBuffer rep ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( rep, getInstanceName() );
		alphabets = XMLParser.extractObjectForTags( xml, "AlphabetContainer", AlphabetContainer.class );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		score = XMLParser.extractObjectForTags( xml, "ScoringFunctions", ScoringFunction[].class );
		set( XMLParser.extractObjectForTags( xml, "index", int[].class ),
				XMLParser.extractObjectForTags( xml, "partialLength", int[].class ),
				XMLParser.extractObjectForTags( xml, "reverse", boolean[].class ));
		try {
			plugIn = XMLParser.extractObjectForTags( xml, "plugIn", boolean.class );
		} catch( Exception e ) {
			plugIn = true;
		}
		extractFurtherInformation( xml );
		setParamsStarts();
	}
	
	/**
	 * This method is used to append further information of the instance to the
	 * XML representation. This method is designed to allow subclasses to add
	 * information to the XML representation.
	 * 
	 * @return the further information as XML code in a {@link StringBuffer}
	 * 
	 * @see AbstractMixtureScoringFunction#extractFurtherInformation(StringBuffer)
	 */
	protected StringBuffer getFurtherInformation() {
		return new StringBuffer( 1 );
	}

	/**
	 * This method is the opposite of {@link #getFurtherInformation()}. It
	 * extracts further information of the instance from a XML representation.
	 * 
	 * @param xml
	 *            the {@link StringBuffer} containing the information to be
	 *            extracted as XML code
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed
	 * 
	 * @see AbstractMixtureScoringFunction#getFurtherInformation()
	 */
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {}


	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getInstanceName()
	 */
	public String getInstanceName() {
		return getClass().getSimpleName();
	}
	
	/**
	 * This method returns a deep copy of the internally used {@link ScoringFunction}.
	 * 
	 * @return a deep copy of the internally used {@link ScoringFunction}
	 * 
	 * @throws Exception if at least one {@link ScoringFunction} could not be cloned
	 * 
	 * @see BasicIndependentProductScoringFunction#getIndices()
	 * @see BasicIndependentProductScoringFunction#getPartialLengths()
	 * @see BasicIndependentProductScoringFunction#getReverseSwitches()
	 */
	public ScoringFunction[] getFunctions() throws Exception {
		return ArrayHandler.clone( score );
	}
	
	/**
	 * This method returns a deep copy of the internally used indices of the {@link ScoringFunction} for the parts.
	 * 
	 * @return a deep copy of the internally used indices of the {@link ScoringFunction} for the parts
	 * 
	 * @see BasicIndependentProductScoringFunction#getFunctions()
	 * @see BasicIndependentProductScoringFunction#getPartialLengths()
	 * @see BasicIndependentProductScoringFunction#getReverseSwitches()
	 */
	public int[] getIndices() {
		return index.clone();
	}
	
	/**
	 * This method returns a deep copy of the internally used partial lengths of the parts.
	 * 
	 * @return a deep copy of the internally used partial lengths of the parts
	 * 
	 * @see BasicIndependentProductScoringFunction#getFunctions()
	 * @see BasicIndependentProductScoringFunction#getIndices()
	 * @see BasicIndependentProductScoringFunction#getReverseSwitches()
	 */
	public int[] getPartialLengths() {
		return partialLength.clone();
	}
	
	/**
	 * This method returns a deep copy of the internally used switches for the parts whether to use the corresponding
	 * {@link ScoringFunction} forward or as reverse complement.
	 * 
	 * @return a deep copy of the internally used switches for the parts whether to use the corresponding 
	 * {@link ScoringFunction} forward or as reverse complement
	 * 
	 * @see BasicIndependentProductScoringFunction#getFunctions()
	 * @see BasicIndependentProductScoringFunction#getIndices()
	 * @see BasicIndependentProductScoringFunction#getPartialLengths()
	 */
	public boolean[] getReverseSwitches() {
		return reverse.clone();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() throws Exception {
		int numPars = this.getNumberOfParameters();
		double[] pars = new double[numPars], help;
		for( int k = 0, i = 0; i < score.length; i++ ) {
			help = score[i].getCurrentParameterValues();
			System.arraycopy( help, 0, pars, k, help.length );
			k += help.length;
		}
		return pars;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScoringFunction#getLogScoreFor(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start ) {
		double s = 0;
		Sequence help;
		for( int myStart, i = 0; i < index.length; i++ ) {
			if( reverse[i] ) {
				try {
					myStart = seq.getLength() - start - this.start[i] - partialLength[i];
					help = seq.reverseComplement();
				} catch ( Exception e ) {
					throw new RuntimeException( e.getMessage() );
				}
			} else {
				help = seq;
				myStart = start + this.start[i];
			}
			
			if( isVariable[index[i]] ) {
				s += ( (VariableLengthScoringFunction)score[index[i]] ).getLogScoreFor( seq, myStart, myStart+partialLength[i]-1 );
			} else {
				s += score[index[i]].getLogScoreFor( seq, myStart );
			}
		}
		return s;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getLogScoreAndPartialDerivation(de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		double s = 0;
		Sequence help;
		for( int myStart, j, i = 0; i < index.length; i++ ) {
			partIList.clear();
			if( reverse[i] ) {
				try {
					myStart = seq.getLength() - start - this.start[i] - partialLength[i];
					help = seq.reverseComplement();
				} catch ( Exception e ) {
					throw new RuntimeException( e.getMessage() );
				}
			} else {
				help = seq;
				myStart = start + this.start[i];
			}
			
			if( isVariable[index[i]] ) {
				s += ( (VariableLengthScoringFunction)score[index[i]] ).getLogScoreAndPartialDerivation( help, myStart, myStart+partialLength[i], partIList, partialDer );
			} else {
				s += score[index[i]].getLogScoreAndPartialDerivation( help, myStart, partIList, partialDer );
			}
			
			for( j = 0; j < partIList.length(); j++ ) {
				indices.add( partIList.get( j ) + params[index[i]] );
			}
		}
		return s;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		if( params == null ) {
			return UNKNOWN;
		} else {
			return params[score.length];
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#getNumberOfRecommendedStarts()
	 */
	public int getNumberOfRecommendedStarts() {
		int max = score[0].getNumberOfRecommendedStarts();
		for( int i = 1; i < score.length; i++ ) {
			max = Math.max( max, score[i].getNumberOfRecommendedStarts() );
		}
		return max;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#setParameters(double[], int)
	 */
	public void setParameters( double[] params, int start ) {
		for( int i = 0; i < score.length; i++ ) {
			score[i].setParameters( params, start + this.params[i] );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 10000 );
		XMLParser.appendObjectWithTags( xml, alphabets, "AlphabetContainer" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, score, "ScoringFunctions" );
		XMLParser.appendObjectWithTags( xml, index, "index" );
		XMLParser.appendObjectWithTags( xml, partialLength, "partialLength" );
		XMLParser.appendObjectWithTags( xml, reverse, "reverse" );
		XMLParser.appendObjectWithTags( xml, plugIn, "plugIn" );
		xml.append( getFurtherInformation() );
		XMLParser.addTags( xml, getInstanceName() );
		return xml;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		StringBuffer sb = new StringBuffer( 100000 );
		for( int i = 0; i < score.length; i++ ) {
			sb.append( "ScoringFunction " + i + ": " + score[i].getInstanceName() + "\n" );
			sb.append( score[i].toString() + "\n" );
		}
		sb.append( "digraph {\n\trankdir=LR\n" );
		sb.append( "\tn" + 0 + "[shape=house, orientation=" + (reverse[0]?90:-90) + ", label=\"emission: " + index[0] + "\\nduration: " + partialLength[0] + "\"]\n" );
		for( int i = 1; i < index.length; i++ ) {
			sb.append( "\tn" + i + "[shape=house, orientation=" + (reverse[i]?90:-90) + ", label=\"emission: " + index[i] + "\\nduration: " + partialLength[i] + "\"]\n" );
			sb.append( "\tn" + (i-1) + "->n" + i  + "\n");
		}
		sb.append( "}" );
		return sb.toString();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#isInitialized()
	 */
	public boolean isInitialized() {
		int i = 0;
		while( i < score.length && score[i].isInitialized() ) {
			i++;
		}
		return i == score.length;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunctionRandomly(boolean)
	 */
	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		for( int i = 0; i < score.length; i++ ) {
			score[i].initializeFunctionRandomly( freeParams );
		}
		setParamsStarts();
	}
}
