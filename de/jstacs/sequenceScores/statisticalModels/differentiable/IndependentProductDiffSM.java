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

package de.jstacs.sequenceScores.statisticalModels.differentiable;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import de.jstacs.NotTrainedException;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.MotifDiscoverer;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.sequenceScores.differentiable.IndependentProductDiffSS;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.AbstractMixtureDiffSM;
import de.jstacs.utils.Normalisation;

/**
 * This class enables the user to model parts of a sequence independent of each
 * other. For instance, the first part of the sequence is modeled by the first
 * {@link DifferentiableStatisticalModel} and has the length of the first
 * {@link DifferentiableStatisticalModel}, the second part starts directly after
 * the first part, is modeled by the second {@link DifferentiableStatisticalModel}
 * ... etc. It is also possible to use a {@link DifferentiableStatisticalModel} for
 * more than one sequence part and in both orientations (if possible).
 * 
 * <br><br>
 * 
 * It is important to set the equivalent sample size (ESS) of each instance carefully, i.e., corresponding to the ESS of the parts.  
 * 
 * @author Jens Keilwagen
 */
public class IndependentProductDiffSM extends IndependentProductDiffSS implements DifferentiableStatisticalModel, MutableMotifDiscoverer
{
	private int[] motifs, componentsUntilNow;
	 /**
	 * The first entry is the index of the corresponding {@link DifferentiableSequenceScore}
	 * in {@link IndependentProductDiffSM#score}. The second entry is
	 * the <b>global</b> motif index in the corresponding
	 * {@link DifferentiableSequenceScore}. The third entry indicates by a &quot;1&quot;
	 * that the corresponding {@link DifferentiableSequenceScore} implements
	 * {@link MutableMotifDiscoverer}, otherwise &quot;-1&quot;.
	 */
	private int[] motifIndexArray;
	/**
	 * The array <code>components</code> stores how many components each {@link DifferentiableSequenceScore} in <code>score</code> has.
	 */
	private int[] components;
	/**
	 * The components of each part.
	 */
	private int[] componentIndexArray;

	private double ess;

	/**
	 * This constructor creates an instance of an
	 * {@link IndependentProductDiffSM} from a given series of
	 * independent {@link DifferentiableStatisticalModel}s. The length that is
	 * modeled by each component is determined by
	 * {@link DifferentiableStatisticalModel#getLength()}. So the length should not be 0.
	 * 
	 * @param ess
	 * 			  the equivalent sample size
	 * @param plugIn whether to use plugIn parameters for the parts, otherwise the last parameters
	 * 			  are used for parts that are instance of
	 * 			  {@link de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM}
	 * @param functions
	 *            the components, i.e. the given series of independent
	 *            {@link DifferentiableStatisticalModel}s
	 *             
	 * @throws CloneNotSupportedException
	 *             if at least one element of <code>functions</code> could not
	 *             be cloned
	 * @throws IllegalArgumentException
	 *             if at least one component has length 0 or if the
	 *             equivalent sample size (ess) is smaller than zero (0)
	 * @throws WrongAlphabetException
	 *             if the user tries to use an alphabet for a reverse complement that can not be used for a reverse complement.
	 *             
	 * @see IndependentProductDiffSM#IndependentProductDiffSM(double, boolean, DifferentiableStatisticalModel[], int[])
	 */
	public IndependentProductDiffSM( double ess, boolean plugIn, DifferentiableStatisticalModel... functions ) throws CloneNotSupportedException,
																						IllegalArgumentException, WrongAlphabetException {
		this( ess, plugIn, functions, IndependentProductDiffSS.getLengthArray( functions ) );
	}

	/**
	 * This constructor creates an instance of an
	 * {@link IndependentProductDiffSM} from given series of
	 * independent {@link DifferentiableStatisticalModel}s and lengths.
	 * 
	 * @param ess
	 * 			  the equivalent sample size
	 * @param plugIn whether to use plugIn parameters for the parts, otherwise the last parameters
	 * 			  are used for parts that are instance of
	 * 			  {@link de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM}
	 * @param functions
	 *            the components, i.e. the given series of independent
	 *            {@link DifferentiableStatisticalModel}s
	 * @param length
	 *            the lengths, one for each component
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one component could not be cloned
	 * @throws IllegalArgumentException
	 *             if the lengths and the components are not matching or if the
	 *             equivalent sample size (ess) is smaller than zero (0)
	 * @throws WrongAlphabetException
	 *             if the user tries to use an alphabet for a reverse complement that can not be used for a reverse complement.
	 * 
	 * @see IndependentProductDiffSM#IndependentProductDiffSM(double, boolean, DifferentiableStatisticalModel[], int[], int[], boolean[])
	 */
	public IndependentProductDiffSM( double ess, boolean plugIn, DifferentiableStatisticalModel[] functions, int[] length ) throws CloneNotSupportedException,
																										IllegalArgumentException, WrongAlphabetException {
		this( ess, plugIn, functions, null, length, null );
	}

	/**
	 * This is the main constructor.
	 *
	 * @param ess the equivalent sample size
	 * @param plugIn whether to use plugIn parameters for the parts, otherwise the last parameters
	 * 			  are used for parts that are instance of
	 * 			  {@link de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM}
	 * @param functions the {@link DifferentiableStatisticalModel}
	 * @param index the index of the {@link DifferentiableStatisticalModel} at each part
	 * @param length the length of each part
	 * @param reverse a switch whether to use it directly or the reverse complementary strand
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one component could not be cloned
	 * @throws IllegalArgumentException
	 *             if the lengths and the components are not matching or if the
	 *             equivalent sample size (ess) is smaller than zero (0)
	 * @throws WrongAlphabetException
	 *             if the user tries to use an alphabet for a reverse complement that can not be used for a reverse complement. 
	 * 
	 */
	public IndependentProductDiffSM( double ess, boolean plugIn, DifferentiableStatisticalModel[] functions, int[] index, int[] length, boolean[] reverse ) throws CloneNotSupportedException,
	IllegalArgumentException, WrongAlphabetException {
		super( plugIn, functions, index, length, reverse );
		this.ess = ess;
		prepare();
	}
	
	private void prepare() {
		motifs = new int[score.length + 1];
		components = new int[score.length];
		for( int i = 0; i < score.length; i++ ) {
			motifs[i + 1] = motifs[i];
			if( score[i] instanceof MotifDiscoverer ) {
				motifs[i + 1] += ( (MotifDiscoverer)score[i] ).getNumberOfMotifs();
			}
			if( score[i] instanceof AbstractMixtureDiffSM ) {
				components[i] = ( (AbstractMixtureDiffSM)score[i] ).getNumberOfComponents();
			} else if( score[i] instanceof MotifDiscoverer ) {
				components[i] = ( (MotifDiscoverer)score[i] ).getNumberOfComponents();
			} else {
				components[i] = 1;
			}
		}
		componentsUntilNow = new int[index.length+1];
		componentsUntilNow[0] = 1;
		for( int i = 0; i < index.length; i++ ){
			componentsUntilNow[i+1] = componentsUntilNow[i] * components[index[i]];
		}
		motifIndexArray = new int[3];
		componentIndexArray = new int[index.length];
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link IndependentProductDiffSM} out of a
	 * {@link StringBuffer} as returned by {@link #toXML()}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public IndependentProductDiffSM( StringBuffer source ) throws NonParsableException {
		super( source );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel#clone()
	 */
	public IndependentProductDiffSM clone() throws CloneNotSupportedException {
		IndependentProductDiffSM clone = (IndependentProductDiffSM)super.clone();
		clone.motifs = motifs.clone();
		clone.motifIndexArray = motifIndexArray.clone();
		clone.components = components.clone();
		clone.componentsUntilNow = componentsUntilNow.clone();
		clone.componentIndexArray = componentIndexArray.clone();
		return clone;
	}

	private int getSFIndex( int index ) {
		int i = 1;
		while( i < startIndexOfParams.length && index >= startIndexOfParams[i] ) {
			i++;
		}
		return i - 1;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		int i = getSFIndex( index );
		return ((DifferentiableStatisticalModel)score[i]).getSizeOfEventSpaceForRandomVariablesOfParameter( index - startIndexOfParams[i] );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogNormalizationConstant()
	 */
	public double getLogNormalizationConstant() {
		double norm = 0;
		for( int i = 0; i < index.length; i++ ) {
			if( isVariable[index[i]] ) {
				norm += ((VariableLengthDiffSM) score[index[i]]).getLogNormalizationConstant( partialLength[i] );
			} else {
				norm += ((DifferentiableStatisticalModel)score[index[i]]).getLogNormalizationConstant();
			}
		}
		return norm;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPartialNormalizationConstant(int)
	 */
	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception {
		int i = 0, j = getSFIndex( parameterIndex ), k = parameterIndex - startIndexOfParams[j];
		double normOfOther = 0, myNorm = 0, myPartNorm = Double.NEGATIVE_INFINITY, n, p;
		for( ; i < index.length; i++ ) {
			if( index[i] == j ) {
				if( isVariable[index[i]] ) {
					p = ((VariableLengthDiffSM) score[j]).getLogPartialNormalizationConstant( k, partialLength[i] );
					n = ((VariableLengthDiffSM) score[j]).getLogNormalizationConstant( partialLength[i] );
				} else {
					p = ((DifferentiableStatisticalModel)score[j]).getLogPartialNormalizationConstant( k );
					n = ((DifferentiableStatisticalModel)score[j]).getLogNormalizationConstant();
				}
				myPartNorm = Normalisation.getLogSum( myPartNorm+n, myNorm+p );
				myNorm += n;
			} else {
				if( isVariable[index[i]] ) {
					normOfOther += ((VariableLengthDiffSM) score[index[i]]).getLogNormalizationConstant( partialLength[i] );
				} else {
					normOfOther += ((DifferentiableStatisticalModel)score[index[i]]).getLogNormalizationConstant();
				}
			}
		}
		return normOfOther + myPartNorm;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getEss()
	 */
	public double getESS() {
		return ess;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.BasicIndependentProductDiffSM#extractFurtherInformation(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInformation( StringBuffer rep ) throws NonParsableException {
		ess = XMLParser.extractObjectForTags( rep, "ess", double.class );
		prepare();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInstanceName()
	 */
	public String getInstanceName() {
		return getClass().getSimpleName();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		if( startIndexOfParams == null ) {
			return UNKNOWN;
		} else {
			return startIndexOfParams[score.length];
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel#getNumberOfRecommendedStarts()
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
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#setParameters(double[], int)
	 */
	public void setParameters( double[] params, int start ) {
		for( int i = 0; i < score.length; i++ ) {
			score[i].setParameters( params, start + this.startIndexOfParams[i] );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.BasicIndependentProductDiffSM#getFurtherInformation()
	 */
	@Override
	protected StringBuffer getFurtherInformation() {
		StringBuffer xml = new StringBuffer( 100 );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
		return xml;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		double val = 0;
		for( int i = 0; i < score.length; i++ ) {
			val += ((DifferentiableStatisticalModel)score[i]).getLogPriorTerm();
		}
		return val;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#addGradientOfLogPriorTerm(double[], int)
	 */
	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		for( int i = 0; i < score.length; i++ ) {
			((DifferentiableStatisticalModel)score[i]).addGradientOfLogPriorTerm( grad, start + startIndexOfParams[i] );
		}
	}	
	
	@Override
	public double getLogProbFor( Sequence sequence, int startpos ) throws Exception {
		return getLogScoreFor( sequence, startpos ) - getLogNormalizationConstant();
	}

	@Override
	public double getLogProbFor( Sequence sequence ) throws Exception {
		return getLogScoreFor( sequence ) - getLogNormalizationConstant();
	}

	@Override
	public double getLogProbFor(Sequence sequence, int startpos, int endpos) throws Exception {
		if( endpos-startpos+1 != length ) {
			throw new IllegalArgumentException();
		} else {
			return getLogProbFor( sequence, startpos );
		}
	}

	/**
	 * This method fills the array {@link IndependentProductDiffSM#motifIndexArray} in the following way.
	 * The first entry is the index of the {@link DifferentiableSequenceScore} that handles the motif with index <code>motifIndex</code>.
	 * The second entry is the index of the motif with index <code>motifIndex</code> in the corresponding {@link DifferentiableSequenceScore}.
	 * The third entry indicates by a &quot;1&quot; that the corresponding {@link DifferentiableSequenceScore} implements
	 * {@link MutableMotifDiscoverer}, otherwise &quot;-1&quot;.
	 * 
	 * @param motifIndex the index of the motif
	 * 
	 * @throws IndexOutOfBoundsException if the motif index is out of bounds
	 * 
	 * @see IndependentProductDiffSM#getNumberOfMotifs()
	 */
	private void fillMotifIndexArray( int motifIndex ) throws IndexOutOfBoundsException {
		motifIndexArray[0] = -1;
		for( int i = 1; i < motifs.length; i++ ) {
			if( score[i - 1] instanceof MotifDiscoverer &&
					( motifIndex >= motifs[i-1] && motifIndex < motifs[i] ) ) {
				motifIndexArray[0] = i - 1;
				motifIndexArray[1] = motifIndex - motifs[motifIndexArray[0]];
				motifIndexArray[2] = score[motifIndexArray[0]] instanceof MutableMotifDiscoverer ? 1 : -1;
				return;
			}
		}
		throw new IndexOutOfBoundsException( "Try to get information about motif " + motifIndex + " (but # motifs = " + getNumberOfMotifs() + ")." );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MutableMotifDiscoverer#initializeMotif(int, de.jstacs.data.DataSet, double[])
	 */
	public void initializeMotif( int motifIndex, DataSet data, double[] weights ) throws Exception {
		fillMotifIndexArray( motifIndex );
		if( motifIndexArray[2] >= 0 ) {
			( (MutableMotifDiscoverer)score[motifIndexArray[0]] ).initializeMotif( motifIndexArray[1], data, weights );
			setParamsStarts();
		} else {
			System.out.println( "warning: not possible" );
		}
	}
	
	public void initializeMotifRandomly( int motif ) throws Exception {
		fillMotifIndexArray( motif );
		if( motifIndexArray[2] >= 0 ) {
			( (MutableMotifDiscoverer)score[motifIndexArray[0]] ).initializeMotifRandomly( motifIndexArray[1] );
			setParamsStarts();
		} else {
			System.out.println( "warning: not possible" );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MutableMotifDiscoverer#modifyMotif(int, int, int)
	 */
	public boolean modifyMotif( int motifIndex, int offsetLeft, int offsetRight ) throws Exception {
		fillMotifIndexArray( motifIndex );
		if( motifIndexArray[2] >= 0 ) {
			boolean b = ( (MutableMotifDiscoverer)score[motifIndexArray[0]] ).modifyMotif( motifIndexArray[1], offsetLeft, offsetRight );
			setParamsStarts();
			return b;
		} else {
			return false;
		}
	}

	/**
	 * This method fills for the user-specified value <code>componentIndex</code> the array
	 * {@link IndependentProductDiffSM#componentIndexArray} with the used component of each
	 * corresponding {@link DifferentiableSequenceScore} of {@link IndependentProductDiffSM#score} using
	 * {@link IndependentProductDiffSM#index}, i.e., it breaks down <code>componentIndex</code>
	 * into the components of each part.
	 * 
	 * @param componentIndex the index of the component in the current model
	 *
	 * @throws IndexOutOfBoundsException if the index is out of bounds
	 */
	private void fillComponentIndexArray( int componentIndex ) throws IndexOutOfBoundsException {
		if( componentIndex < 0 || componentIndex > getNumberOfComponents() ) {
			throw new IndexOutOfBoundsException();
		}
		for( int i = index.length - 1; i >= 0; i-- ) {
			componentIndexArray[i] = componentIndex % components[i];
			componentIndex /= components[i];
		}
	}

	/**
	 * This method fills for the user-specified values <code>componentIndex</code> and <code>motifIndex</code>
	 * the arrays {@link IndependentProductDiffSM#motifIndexArray} and {@link IndependentProductDiffSM#componentIndexArray}.
	 * For the later see method {@link IndependentProductDiffSM#fillComponentIndexArray(int)}.
	 * For {@link IndependentProductDiffSM#motifIndexArray}, it sets the with the following values. 
	 * The first entry is the index of the corresponding {@link DifferentiableSequenceScore} in {@link IndependentProductDiffSM#score}.
	 * The second entry is the <b>global</b> motif index in the corresponding {@link DifferentiableSequenceScore}.
	 * The third entry indicates by a &quot;1&quot; that the corresponding {@link DifferentiableSequenceScore} implements
	 * {@link MutableMotifDiscoverer}, otherwise &quot;-1&quot;.
	 * 
	 * @param componentIndex the index of the component of the current instance
	 * @param motifIndex the (local) index of the motif of the current instance
	 * 
	 * @throws IndexOutOfBoundsException if the indices are out of bounds
	 */
	private void fillMotifIndexArrayForComponent( int componentIndex, int motifIndex ) throws IndexOutOfBoundsException {
		fillComponentIndexArray( componentIndex );
		int i = -1, anz = 0;
		do {
			motifIndex -= anz;
			i++;
			if( i == score.length ) {
				throw new IndexOutOfBoundsException();
			}
			fillHashWith( i );
			anz = hash.size();
		} while( anz != 0 && motifIndex < anz );
		motifIndexArray[0] = i;
		motifIndexArray[2] = score[motifIndexArray[0]] instanceof MutableMotifDiscoverer ? 1 : -1;
		
		int[] specificMotifs = new int[anz];
		Iterator<Integer> it = hash.iterator();
		for( i = 0; i < anz; i++ ) {
			specificMotifs[i] = it.next();
		}
		Arrays.sort( specificMotifs );
		
		motifIndexArray[1] = specificMotifs[motifIndex];
	}
	
	private HashSet<Integer> hash;
	
	/**
	 * This method collects in the {@link HashSet} {@link IndependentProductDiffSM#hash}
	 * the global indices of the motifs of {@link DifferentiableSequenceScore} <code>score[index]</code> for a given filled
	 * array {@link IndependentProductDiffSM#components}.
	 *  
	 * @param scoreIndex the index of the {@link DifferentiableSequenceScore}.
	 * 
	 * @see IndependentProductDiffSM#score
	 */
	private void fillHashWith( int scoreIndex ) {
		if ( hash == null ) {
			hash = new HashSet<Integer>();
		}
		if( scoreIndex == score.length ) {
			throw new IndexOutOfBoundsException();
		}
		hash.clear();
		if( score[scoreIndex] instanceof MotifDiscoverer ) {
			MotifDiscoverer md = (MotifDiscoverer) score[scoreIndex];
			for( int h, l, m, idx = 0; idx < index.length; idx++ ) {
				if( index[idx] == scoreIndex ) {
					m = md.getNumberOfMotifsInComponent( components[idx] );
					for( l = 0; l < m; l++ ) {
						h = md.getGlobalIndexOfMotifInComponent( components[idx], l );
						if( !hash.contains( h ) ) {
							hash.add( h );
						}
					}
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getGlobalIndexOfMotifInComponent(int, int)
	 */
	public int getGlobalIndexOfMotifInComponent( int component, int motif ) {
		fillMotifIndexArrayForComponent( component, motif );
		return motifs[motifIndexArray[0]] + motifIndexArray[1];
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getIndexOfMaximalComponentFor(de.jstacs.data.Sequence)
	 */
	public int getIndexOfMaximalComponentFor( Sequence sequence ) throws Exception {
		int c = 0, idx;
		for( int i = 0; i < index.length; i++ ) {
			if( score[index[i]] instanceof AbstractMixtureDiffSM ) {
				if( reverse[i] ) {
					idx = ( (AbstractMixtureDiffSM)score[index[i]] ).getIndexOfMaximalComponentFor(
							sequence.reverseComplement(), sequence.getLength() - start[i] - partialLength[i] );
				} else {
					idx = ( (AbstractMixtureDiffSM)score[index[i]] ).getIndexOfMaximalComponentFor( sequence, start[i] );
				}
			} else if( score[i] instanceof MotifDiscoverer ) {
				if( reverse[i] ) {
					idx = ( (MotifDiscoverer)score[index[i]] ).getIndexOfMaximalComponentFor(
							sequence.reverseComplement( sequence.getLength() - start[i] - partialLength[i],
									sequence.getLength() - start[i] ) );
				} else {
					idx = ( (MotifDiscoverer)score[index[i]] ).getIndexOfMaximalComponentFor(
							sequence.getSubSequence( score[i].getAlphabetContainer(),
							start[i],
							partialLength[i] ) );
				}
			} else {
				idx = 0;
			}
			c += componentsUntilNow[i] * idx;
		}
		return c;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getMotifLength(int)
	 */
	public int getMotifLength( int motif ) {
		fillMotifIndexArray( motif );
		return ( (MotifDiscoverer)score[motifIndexArray[0]] ).getMotifLength( motifIndexArray[1] );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getNumberOfComponents()
	 */
	public int getNumberOfComponents() {
		return componentsUntilNow[index.length];
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getNumberOfMotifs()
	 */
	public int getNumberOfMotifs() {
		return motifs[score.length];
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getNumberOfMotifsInComponent(int)
	 */
	public int getNumberOfMotifsInComponent( int component ) {
		fillComponentIndexArray( component );
		int i = 0, anz = 0;
		for( ; i < score.length; i++ ) {
			fillHashWith( i );
			anz += hash.size();
		}
		return anz;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getProfileOfScoresFor(int, int, de.jstacs.data.Sequence, int, de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile)
	 */
	public double[] getProfileOfScoresFor( int component, int motif, Sequence sequence, int startpos, KindOfProfile dist ) throws Exception {
		fillMotifIndexArrayForComponent( component, motif );
		MotifDiscoverer md = (MotifDiscoverer)score[motifIndexArray[0]];
		int w =  md.getMotifLength( motifIndexArray[1] );
		double[] erg = new double[length - w + 1];
		double[] part;
		Arrays.fill( erg, Double.NEGATIVE_INFINITY );
		for( int idx, k, j, i = 0; i < index.length; i++ ) {
			if( index[i] == motifIndexArray[0] ) {
				// part from that DifferentiableSequenceScore that contains the motif
				
				//find local index
				idx = 0;
				while( md.getGlobalIndexOfMotifInComponent( componentIndexArray[i], idx ) != motifIndexArray[1] ) {
					idx++;
				}
				//compute profile and fill in result array erg
				if( reverse[i] ) {
					part = md.getProfileOfScoresFor( componentIndexArray[i], idx,
							sequence.reverseComplement(), sequence.getLength() - startpos - start[i] - partialLength[i], dist );
					System.arraycopy( part, 0, erg, start[i], partialLength[i] - w + 1 );
				} else {
					part = md.getProfileOfScoresFor( componentIndexArray[i], idx,
							sequence, startpos + start[i], dist );
					for( k = start[i], j = part.length-1; j >= 0; j--, k++ ){
						erg[k] = part[j];
					}
				}
			}
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getStrandFor(int, int, de.jstacs.data.Sequence, int)
	 */
	public double[] getStrandProbabilitiesFor( int component, int motif, Sequence sequence, int startpos ) throws Exception {
		//TODO
		return new double[]{0.5,0.5};
	}
	
	
	public boolean isNormalized() {
		return AbstractDifferentiableStatisticalModel.isNormalized( score );
	}

	public void adjustHiddenParameters( int index, DataSet[] data, double[][] weights ) throws Exception {
		DataSet[] part = new DataSet[data.length];
		double[][] help;
		for( int a, i = 0; i < score.length; i++ ) {
			if( score[i] instanceof MutableMotifDiscoverer ) {
				a = extractSequenceParts( i, data, part );
				help = a==1 ? weights : extractWeights( a, weights );
				((MutableMotifDiscoverer) score[i]).adjustHiddenParameters( index, part, help );
			}
		}		
		setParamsStarts();
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.StatisticalModel#emitDataSet(int, int[])
	 */
	@Override
	public DataSet emitDataSet(int numberOfSequences, int... seqLength) throws NotTrainedException, Exception {
		throw new Exception( "Standard implementation of emitDataSet used for "
						+ getInstanceName()	+ ". You have to overwrite this method to use it in a proper way.");
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.StatisticalModel#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() throws UnsupportedOperationException {
		throw new UnsupportedOperationException( "The maximal markov order for this model in undefined.");
	}
}
