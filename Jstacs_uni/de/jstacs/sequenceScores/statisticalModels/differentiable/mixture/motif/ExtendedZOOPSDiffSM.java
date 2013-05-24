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

package de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.motif;

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.Mutable;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.NormalizedDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.AbstractMixtureDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;

/**
 * This class handles mixtures with at least one hidden motif.
 * 
 * @author Jens Keilwagen
 */
public class ExtendedZOOPSDiffSM extends AbstractMixtureDiffSM implements MutableMotifDiscoverer
{
	/**
	 * This constant indicates that in each sequence has one binding site of a motif instance (similar to OOPS).
	 */
	public static final boolean CONTAINS_ALWAYS_A_MOTIF = true;

	/**
	 * This constant indicates that a sequence possibly has one binding site of a motif instance (similar to ZOOPS).
	 */
	public static final boolean CONTAINS_SOMETIMES_A_MOTIF = false;

	private double[][] simpleScore;

	private double[] bgHelp;

	private int[] anz, currentPos;

	private double ess;

	private int bgIndex;

	private boolean type;
	
	/**
	 * Indicates whether the old parameters or new one should be used.
	 */
	private boolean plugInBg;

	private HomogeneousDiffSM bg;

	private static DifferentiableStatisticalModel[] getDifferentiableStatisticalModels( int length,
			HomogeneousDiffSM bg, DifferentiableStatisticalModel[] motif, DurationDiffSM[] posPrior )
	{
		int m = motif.length;
		if( m == 0 )
		{
			throw new IllegalArgumentException(	"Please insert at least one DifferentiableSequenceScore for a motif." );
		}
		DifferentiableStatisticalModel[] erg = new DifferentiableStatisticalModel[2 * m + 1];
		if( posPrior == null )
		{
			posPrior = new DurationDiffSM[m];
		} else if ( m != posPrior.length ) {
			throw new IllegalArgumentException(	"The number of motifs and durations has to be the same." );
		}
		int i = 0, l = -1;
		AlphabetContainer seqABC = bg.getAlphabetContainer();
		if( !seqABC.isSimple() )
		{
			throw new IllegalArgumentException(
					"The AlphabetContainer of the motif and background models has to be simple." );
		}
		AlphabetContainer posABC = null;
		while( i < m )
		{
			if( !seqABC.checkConsistency( motif[i].getAlphabetContainer() ) )
			{
				throw new IllegalArgumentException( "The " + i
						+ "-th motif model is not correct. Check the AlphabetContainer." );
			}
			erg[2 * i] = motif[i];
			if( posPrior[i] != null )
			{
				if( 0 > posPrior[i].getMin() || posPrior[i].getMax() > length - motif[i].getLength() )
				{
					throw new IllegalArgumentException( "The " + i
							+ "-th model for the positional information is not correct. Check the AlphabetContainer." );
				}
				erg[2 * i + 1] = posPrior[i];
			}
			else
			{
				if( posABC == null || motif[i].getLength() != l )
				{
					posABC = new AlphabetContainer( new DiscreteAlphabet( 0, length - motif[i].getLength() ) );
				}
				erg[2 * i + 1] = new UniformDurationDiffSM( 0, length - motif[i].getLength() );
				l = motif[i].getLength();
			}
			i++;
		}
		if( !seqABC.checkConsistency( bg.getAlphabetContainer() ) )
		{
			throw new IllegalArgumentException( "The background model is not correct. Check the AlphabetContainer." );
		}
		erg[2 * i] = bg;
		return erg;
	}
	

	/**
	 * This constructor creates an instance of {@link ExtendedZOOPSDiffSM} that is either an OOPS or a ZOOPS model depending on the chosen <code>type</code>.
	 * 
	 * @param type the type of hidden motifs, either {@link ExtendedZOOPSDiffSM#CONTAINS_ALWAYS_A_MOTIF} or {@link ExtendedZOOPSDiffSM#CONTAINS_SOMETIMES_A_MOTIF}
	 * @param length the length of the modeled sequences (e.g. 500)
	 * @param starts the number of recommended starts
	 * @param plugIn a switch whether to use plug-in or randomly chosen parameter when using {@link DifferentiableSequenceScore#initializeFunction(int, boolean, DataSet[], double[][])}
	 * @param bg the {@link DifferentiableSequenceScore} for the overall background (i.e. flanking sequence that does not contain a motif)
	 * @param motif the {@link DifferentiableSequenceScore} for the motif
	 * @param posPrior the {@link DifferentiableSequenceScore} for the position
	 * @param plugInBg a switch whether to plug in the (old = last) parameters of the background model, if <code>false</code> the background model is initialized again
	 * 
	 * @throws Exception if something went wrong (e.g. the {@link AlphabetContainer} is not simple, ...)
	 */
	public ExtendedZOOPSDiffSM( boolean type, int length, int starts, boolean plugIn, HomogeneousDiffSM bg,
			DifferentiableStatisticalModel motif, DurationDiffSM posPrior, boolean plugInBg ) throws Exception
	{
		this( type, length, starts, plugIn, bg, new DifferentiableStatisticalModel[]{ motif }, new DurationDiffSM[]{ posPrior }, plugInBg );
	}
	
	/**
	 * This constructor creates an instance of {@link ExtendedZOOPSDiffSM} that allows to have one site of the specified motifs in a {@link Sequence}.
	 * 
	 * @param type the type of hidden motifs, either {@link ExtendedZOOPSDiffSM#CONTAINS_ALWAYS_A_MOTIF} or {@link ExtendedZOOPSDiffSM#CONTAINS_SOMETIMES_A_MOTIF}
	 * @param length the length of the modeled sequences (e.g. 500)
	 * @param starts the number of recommended starts
	 * @param plugIn a switch whether to use plug-in or randomly chosen parameter when using {@link DifferentiableSequenceScore#initializeFunction(int, boolean, DataSet[], double[][])}
	 * @param bg the {@link DifferentiableSequenceScore} for the overall background (i.e. flanking sequence that does not contain a motif)
	 * @param motif the {@link DifferentiableSequenceScore}s for the sequence motif
	 * @param posPrior the {@link DifferentiableSequenceScore}s for the position of the the sequence motifs
	 * @param plugInBg a switch whether to plug in the parameters of the background model in some way, if <code>false</code> the background model is initialized uniformly
	 * 
	 * @throws Exception if something went wrong (e.g. the {@link AlphabetContainer} is not simple, ...)
	 */
	public ExtendedZOOPSDiffSM( boolean type, int length, int starts, boolean plugIn, HomogeneousDiffSM bg,
			DifferentiableStatisticalModel[] motif, DurationDiffSM[] posPrior, boolean plugInBg )
			throws Exception
	{
		super( length, starts, motif.length + (type == CONTAINS_ALWAYS_A_MOTIF ? 0 : 1), true, plugIn,
				getDifferentiableStatisticalModels( length, bg, motif, posPrior ) );
		if( !alphabets.isSimple() )
		{
			throw new IllegalArgumentException( "The AlphabetContainer has to be simple." );
		}
		this.type = type;
		this.plugInBg = plugInBg;
		initObject();
		setHyperParametersOfBackgroundModel();
		computeLogGammaSum();
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates a new {@link ExtendedZOOPSDiffSM} out of a
	 * {@link StringBuffer} as returned by {@link #toXML()}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the representation could not be parsed
	 */
	public ExtendedZOOPSDiffSM( StringBuffer source ) throws NonParsableException
	{
		super( source );
		initObject();
	}

	private void initObject()
	{
		bgIndex = function.length - 1;
		bg = (HomogeneousDiffSM) function[bgIndex];

		simpleScore = new double[(function.length-1)/2][];
		if( type == CONTAINS_ALWAYS_A_MOTIF )
		{
			ess = 0;
		}
		else
		{
			ess = bg.getESS();
		}
		for( int i = 0; i < bgIndex; i += 2 )
		{
			ess += function[i].getESS();
		}
		createSimpleScore();
		initBgHelp();
		anz = new int[componentScore.length + 1];
		currentPos = new int[1];
	}

	private void createSimpleScore()
	{
		for( int i = 0; i < bgIndex; i += 2 )
		{
			simpleScore[i / 2] = new double[(int) function[i + 1].getAlphabetContainer().getAlphabetLengthAt( 0 )];
		}
	}
	
	private void setHyperParametersOfBackgroundModel() throws Exception
	{
		int[] len = new int[length+1], l = new int[1];
		int i = 0, anz, m;
		for( ; i <= length; i++ )
		{
			len[i] = i;
		}
		double[] weight = new double[len.length];
		double w;
		boolean okay;
		DurationDiffSM dur;
		for( i = 0; i < bgIndex; i += 2 )
		{
			m = length - function[i].getLength();
			dur = (DurationDiffSM) function[i+1];
			
			//determine how many possibilities can be used for the current motif
			anz = 0;
			dur.reset();
			do
			{
				dur.getInternalPosition(l);
				if( l[0] <= m )
				{
					anz++;
					okay = dur.next();
				}
				else
				{
					okay = false;
				}
			}
			while( okay );
			
			//determine the weight for each possible length
			w = function[i].getESS() / (double) anz;
			
			//create the statistic
			dur.reset();
			do
			{
				dur.getInternalPosition(l);
				if( l[0] <= m )
				{
					// add the weight for the sequence before and after the motif
					weight[l[0]] += w;
					weight[m - l[0]] += w;
					okay = dur.next();
				}
				else
				{
					okay = false;
				}
			}
			while( okay );			
		}
		
		// add the weight if no motif is in the sequence
		if( type == CONTAINS_SOMETIMES_A_MOTIF )
		{
			weight[length] += bg.getESS();
		}
		
		bg.setStatisticForHyperparameters( len, weight );
	}
	
	public ExtendedZOOPSDiffSM clone() throws CloneNotSupportedException
	{
		ExtendedZOOPSDiffSM clone = (ExtendedZOOPSDiffSM) super.clone();
		clone.simpleScore = new double[simpleScore.length][];
		for( int i = 0; i < simpleScore.length; i++ )
		{
			clone.simpleScore[i] = new double[simpleScore[i].length];
		}
		clone.bg = (HomogeneousDiffSM) clone.function[bgIndex];
		clone.bgHelp = bgHelp.clone();
		clone.anz = anz.clone();
		clone.currentPos = currentPos.clone();
		return clone;
	}

	public double getHyperparameterForHiddenParameter( int index )
	{
		return function[2 * index].getESS();
	}

	protected double getLogNormalizationConstantForComponent( int i )
	{
		i *= 2;
		if( i >= 0 && i < bgIndex )
		{
			int l = length - function[i].getLength();
			//System.out.println( function[i].getLogNormalizationConstant() + " + " + function[i + 1].getLogNormalizationConstant() + " + " + bg.getLogNormalizationConstant( l ) );
			return function[i].getLogNormalizationConstant()
				+ function[i + 1].getLogNormalizationConstant()
				+ bg.getLogNormalizationConstant( l );
		}
		else if( i == bgIndex )
		{
			return bg.getLogNormalizationConstant();
		}
		else
		{
			throw new IndexOutOfBoundsException();
		}

	}

	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception
	{
		if( isNormalized() )
		{
			return Double.NEGATIVE_INFINITY;
		}
		else
		{
			if( Double.isNaN( norm ) )
			{
				precomputeNorm();
			}
			int[] ind = getIndices( parameterIndex );
			double res;
			if( ind[0] == function.length )
			{
				res = partNorm[ind[1]];
			}
			else if( ind[0] == bgIndex )
			{
				res = Double.NEGATIVE_INFINITY;
				for( int l, i = 0; i < bgIndex; i += 2 )
				{
					l = length - function[i].getLength();
					res = Normalisation.getLogSum( res, logHiddenPotential[i / 2] + function[i].getLogNormalizationConstant()
							+ function[i + 1].getLogNormalizationConstant() + bg.getLogPartialNormalizationConstant( ind[1], l ) );
				}
				if( type == CONTAINS_SOMETIMES_A_MOTIF )
				{
					res = Normalisation.getLogSum( res, logHiddenPotential[bgIndex / 2] + bg.getLogPartialNormalizationConstant( ind[1] ) );
				}
			}
			else
			{
				if( ind[0] % 2 == 0 )
				{
					int l = length - function[ind[0]].getLength();
					res = logHiddenPotential[ind[0] / 2]
					      + function[ind[0]].getLogPartialNormalizationConstant( ind[1] )
					      + function[ind[0] + 1].getLogNormalizationConstant()
					      + bg.getLogNormalizationConstant( l );
				}
				else
				{
					int l = length - function[ind[0] - 1].getLength();
					res = logHiddenPotential[ind[0] / 2]
					      + function[ind[0] - 1].getLogNormalizationConstant()
					      + function[ind[0]].getLogPartialNormalizationConstant( ind[1] )
					      + bg.getLogNormalizationConstant( l );
				}
			}
			
			//System.out.println( parameterIndex + "\t" + Arrays.toString(ind) + "\t" + res );
			
			return res;
		}
	}

	public double getESS()
	{
		return ess;
	}

	public void initializeFunctionRandomly( boolean freeParams ) throws Exception
	{
		for( int n = function.length-1, i = 0; i < n; i++ )
		{
			function[i].initializeFunctionRandomly( freeParams );
		}
		if( !plugInBg ) {
			bg.initializeFunctionRandomly( freeParams );
		}
		if( optimizeHidden ) {
			initializeHiddenPotentialRandomly();
		}
		init( freeParams );
		initBgHelp();
	}
	
	private void initBgHelp()
	{
		int n = bg.getNumberOfParameters();
		if( n == DifferentiableSequenceScore.UNKNOWN )
		{
			bgHelp = new double[0];
		}
		else if( bgHelp == null || bgHelp.length != n )
		{
			bgHelp = new double[n];
		}
	}
	
	public String getInstanceName()
	{
		String erg = "hiddenMotifsMixture("
				+ (type == CONTAINS_ALWAYS_A_MOTIF ? "contains always a motif" : "contains sometimes a motif") + "; bg = "
				+ bg.getInstanceName();
		for( int i = 0; i < function.length-1; i += 2 )
		{
			erg += "; (" + function[i].getInstanceName() + ", " + function[i+1].getInstanceName() + ")";
		}
		return erg + ")";
	}

	/**
	 * This method fills an internal array with the partial scores.
	 * 
	 * @param i the index of the component
	 * @param seq the {@link Sequence}
	 * @param start the start position
	 * 
	 * @return the number of computed partial scores
	 */
	protected int fillComponentScoreOf( int i, Sequence seq, int start )
	{
		int j = 2 * i, m = function[j].getLength(), l = 0, bgOrder = bg.getMaximalMarkovOrder();
		//start & end points
		int homSt, homE;
		PositionDiffSM pos = (PositionDiffSM) function[j + 1];
		pos.reset();
		do
		{
			pos.getInternalPosition( currentPos );
			homSt = Math.max( 0, currentPos[0] - bgOrder );
			homE = Math.min( length, currentPos[0] + m + bgOrder )-1;
			simpleScore[i][l++] = pos.getLogScoreForInternal() + function[j].getLogScoreFor( seq, start + currentPos[0] )
				- bg.getLogScoreFor( seq, start + homSt, homE )
				+ bg.getLogScoreFor( seq, start + homSt, currentPos[0]-1 ) // left
				+ bg.getLogScoreFor( seq, start + currentPos[0] + m, homE ); //right

		} while( pos.next() );
		return l;
	}

	protected void fillComponentScores( Sequence seq, int start )
	{
		int i = 0, j;
		for( ; i < bgIndex / 2; i++ )
		{
			j = fillComponentScoreOf( i, seq, start );
			componentScore[i] = logHiddenPotential[i] + Normalisation.getLogSum( 0, j, simpleScore[i] );
		}
		if( type == CONTAINS_SOMETIMES_A_MOTIF )
		{
			componentScore[i] = logHiddenPotential[i];
		}
	}

	public double getLogScoreFor( Sequence seq, int start )
	{
		fillComponentScores( seq, start );
		return bg.getLogScoreFor( seq, start, start+length-1 ) + Normalisation.getLogSum( componentScore );
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer )
	{
		int i = 0, j = 0, l, m, n, counter, stop, bgOrder = bg.getMaximalMarkovOrder();
		int homSt, homE;
		anz[0] = partialDer.length();
		int[][] end = new int[4][];
		PositionDiffSM pos;
		// for each motif (function)
		for( ; j < bgIndex; i++, j = 2 * i )
		{
			m = function[j].getLength();
			stop = length - m + 1;
			
			iList[j].clear();
			dList[j].clear();
			end[0] = new int[stop];
			iList[j + 1].clear();
			dList[j + 1].clear();
			pos = (PositionDiffSM) function[j + 1];
			pos.reset();
			end[1] = new int[stop];
			iList[bgIndex].clear();
			dList[bgIndex].clear();
			end[2] = new int[stop];
			end[3] = new int[stop];
			
			// for each start position
			stop = 0;
			do
			{
				// get current position
				pos.getInternalPosition( currentPos );
				homSt = Math.max( 0, currentPos[0] - bgOrder );
				homE = Math.min( length, currentPos[0] + m + bgOrder )-1;
				// compute the score
				simpleScore[i][stop] =
					pos.getLogScoreAndPartialDerivationForInternal( iList[j + 1], dList[j + 1] ) // position
					+ function[j].getLogScoreAndPartialDerivation( seq, start + currentPos[0], iList[j], dList[j] ) // motif
					- bg.getLogScoreAndPartialDerivation( seq, start + homSt, homE, iList[bgIndex], dList[bgIndex] );
				
				end[0][stop] = iList[j].length();
				end[1][stop] = iList[j + 1].length();
				end[2][stop] = iList[bgIndex].length();
				
				simpleScore[i][stop] +=
					bg.getLogScoreAndPartialDerivation( seq, start + homSt, currentPos[0]-1, iList[bgIndex], dList[bgIndex] )
					+ bg.getLogScoreAndPartialDerivation( seq, start + currentPos[0] + m, homE, iList[bgIndex], dList[bgIndex] );

				end[3][stop++] = iList[bgIndex].length();
			}while( pos.next() );

			// normalize
			componentScore[i] = logHiddenPotential[i]
					+ Normalisation.logSumNormalisation( simpleScore[i], 0, stop, simpleScore[i], 0 );

			// weight the partial derivations
			for( m = 0; m < 2; m++ )
			{
				for( l = 0, counter = 0; l < stop; l++ )
				{
					while( counter < end[m][l] )
					{
						indices.add( iList[j + m].get( counter ) + paramRef[j + m] );
						partialDer.add( dList[j + m].get( counter++ ) * simpleScore[i][l] );
					}
				}
			}
			Arrays.fill( bgHelp, 0 );
			for( l = 0, counter = 0; l < stop; l++ )
			{
				while( counter < end[2][l] )
				{
					bgHelp[iList[bgIndex].get( counter )] -= (dList[bgIndex].get( counter ) * simpleScore[i][l]);
					counter++;
				}
				while( counter < end[3][l] )
				{
					bgHelp[iList[bgIndex].get( counter )] += (dList[bgIndex].get( counter ) * simpleScore[i][l]);
					counter++;
				}
			}
			for( counter = 0; counter < bgHelp.length; counter++ )
			{
				indices.add( paramRef[bgIndex] + counter );
				partialDer.add( bgHelp[counter] );
			}
			anz[i + 1] = partialDer.length();
		}
		if( type == CONTAINS_SOMETIMES_A_MOTIF )
		{
			componentScore[bgIndex / 2] = logHiddenPotential[bgIndex / 2];
		}

		double logScore = Normalisation.logSumNormalisation( componentScore, 0, componentScore.length, componentScore, 0 );

		// adjust if necessary
		if( componentScore.length > 1 )
		{
			for( i = 0; i < hiddenPotential.length; i++ )
			{
				partialDer.multiply( anz[i], anz[i + 1], componentScore[i] );
			}
		}

		// bg for complete sequence
		iList[bgIndex].clear();
		dList[bgIndex].clear();
		logScore += bg.getLogScoreAndPartialDerivation( seq, start, start+length-1, iList[bgIndex], dList[bgIndex] );

		for( i = 0; i < iList[bgIndex].length(); i++ )
		{
			indices.add( paramRef[bgIndex] + iList[bgIndex].get( i ) );
			partialDer.add( dList[bgIndex].get( i ) );
		}

		// hiddenLambda
		n = bgIndex + 1;
		i = paramRef[n + 1] - paramRef[n];
		for( j = 0; j < i; j++ )
		{
			indices.add( paramRef[n] + j );
			partialDer.add( componentScore[j] - (isNormalized()?hiddenPotential[j]:0) );
		}

		return logScore;
	}

	public String toString( NumberFormat nf )
	{
		if( Double.isNaN( norm ) )
		{
			precomputeNorm();
		}
		StringBuffer erg = new StringBuffer( function.length * 1000 );
		erg.append( "bg:\n" + bg.toString(nf) + "\n" );
		if( type == CONTAINS_SOMETIMES_A_MOTIF )
		{
			erg.append( "\nno motif: " + nf.format(Math.exp(partNorm[bgIndex / 2]) - norm) + "\texp(" +partNorm[bgIndex / 2] + " - " + norm + ")\t" + logHiddenPotential[bgIndex/2] + "\n" );
		}
		for( int i = 0; i < bgIndex; i += 2 )
		{
			erg.append( "\nmotif " + (i / 2) + ": " );
			if( hiddenPotential.length > 1 )
			{
				nf.format(erg.append( Math.exp(partNorm[i / 2] - norm)) +"\texp(" +partNorm[i / 2] + " - " + norm + ")\t" + logHiddenPotential[i/2] );
			}
			erg.append( "\n" + function[i].toString(nf) + "\n" + function[i + 1].toString(nf) + "\n" );
		}
		return erg.toString();
	}

	protected StringBuffer getFurtherInformation()
	{
		StringBuffer erg = new StringBuffer( 200 );
		XMLParser.appendObjectWithTags( erg, type, "type" );
		XMLParser.appendObjectWithTags( erg, plugInBg, "plugInBg" );
		return erg;
	}

	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException
	{
		type = XMLParser.extractObjectForTags( xml, "type", boolean.class );
		plugInBg = XMLParser.extractObjectForTags( xml, "plugInBg", boolean.class );
	}
	
	public boolean modifyMotif( int motif, int offsetLeft, int offsetRight ) throws Exception
	{		
		int c = getNumberOfComponents() - 1;
		if( (motif < c || (motif == c && type == CONTAINS_ALWAYS_A_MOTIF) )
				&& function[2*motif] instanceof Mutable )
		{
			double norm_old = function[2*motif].getLogNormalizationConstant();
			boolean res = ((Mutable) function[2*motif]).modify(offsetLeft, offsetRight);
			if( res )
			{
				setHyperParametersOfBackgroundModel();
				init( freeParams );
				((DurationDiffSM)function[2*motif+1]).modify( offsetLeft-offsetRight );
				simpleScore[motif] = new double[(int) function[2*motif + 1].getAlphabetContainer().getAlphabetLengthAt( 0 )];
				double norm_new = function[2*motif].getLogNormalizationConstant();
				hiddenParameter[motif] += ( norm_old - norm_new );
				norm = Double.NaN;
				this.setHiddenParameters( hiddenParameter , 0 );
			}
			
			//System.out.println( getNormalizationConstant() + "\t" + getNormalizationConstantForComponent( motif ) + "\t" + hiddenParameter[motif] );
			return res;
		} else {
			return false;
		}
	}

	public void initializeMotif( int motif, DataSet data, double[] weights ) throws Exception
	{
		int c = getNumberOfComponents() - 1;
		if( motif < c || (motif == c && type == CONTAINS_ALWAYS_A_MOTIF) )
		{
			function[2*motif].initializeFunction( 0, freeParams, new DataSet[]{data}, (weights==null ? null : new double[][]{weights}) );
			init( freeParams );
		}
		else
		{
			throw new IndexOutOfBoundsException();
		}
	}
	
	public void initializeMotifRandomly( int motif ) throws Exception
	{
		int c = getNumberOfComponents() - 1;
		if( motif < c || (motif == c && type == CONTAINS_ALWAYS_A_MOTIF) )
		{
			function[2*motif].initializeFunctionRandomly(freeParams);
			function[2*motif+1].initializeFunctionRandomly(freeParams);
			init( freeParams );
		}
		else
		{
			throw new IndexOutOfBoundsException();
		}
	}

	public int getNumberOfMotifs()
	{
		return getNumberOfComponents() - (type==CONTAINS_ALWAYS_A_MOTIF?0:1);
	}

	public int getNumberOfMotifsInComponent( int component )
	{
		int c = getNumberOfComponents() - 1;
		if( component < c || (component == c && type == CONTAINS_ALWAYS_A_MOTIF) )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	
	public int getIndexOfMaximalComponentFor( Sequence sequence )
	{
		return getIndexOfMaximalComponentFor( sequence, 0 );
	}

	public int getGlobalIndexOfMotifInComponent( int component, int motif )
	{
		return component;
	}

	public double[] getProfileOfScoresFor( int component, int motif, Sequence sequence, int startpos, KindOfProfile dist ) throws WrongLengthException
	{
		if( length != sequence.getLength() - startpos )
		{
			throw new WrongLengthException( "The model can not score a sequence of this length." );
		}
		int c = getNumberOfComponents() - (type==CONTAINS_SOMETIMES_A_MOTIF?1:0);
		if( motif == 0 && component < c )
		{
			//unnormalized log score
			fillComponentScoreOf( component, sequence, startpos );

			double d = 0;
			switch( dist )
			{
				case UNNORMALIZED_JOINT:
					d = logHiddenPotential[component];
				case UNNORMALIZED_CONDITIONAL:
					d += bg.getLogScoreFor( sequence, startpos, startpos+length-1 );
					break;				
				case NORMALIZED_CONDITIONAL:
					d = -Normalisation.getLogSum( 0, simpleScore[component].length, simpleScore[component] );
					break;
				
				default:
					throw new IndexOutOfBoundsException();
			}
			
			double[] res = new double[length - function[2*component].getLength() + 1];
			DurationDiffSM posPrior = (DurationDiffSM) function[2*component+1];
			posPrior.reset();
			int[] internal = new int[1];
			posPrior.getInternalPosition( internal );
			boolean b;
			for( int i = 0, j = 0; i < res.length; i++ )
			{
				if( i == internal[0] )
				{
					res[i] = simpleScore[component][j++] + d;
					b = posPrior.next();
					if( b )
					{
						posPrior.getInternalPosition(internal);
					}
					else
					{
						internal[0] = -1;
					}
				}
				else
				{
					res[i] = Double.NEGATIVE_INFINITY;
				}
			}
			return res;
		}
		else
		{
			throw new IndexOutOfBoundsException();
		}
	}

	public int getMotifLength( int motif )
	{
		int c = getNumberOfComponents() - 1;
		if( motif < c || (motif == c && type == CONTAINS_ALWAYS_A_MOTIF) )
		{
			return function[2*motif].getLength();
		}
		else
		{
			throw new IndexOutOfBoundsException();
		}
	}	

	public void adjustHiddenParameters( int classIndex, DataSet[] data, double[][] dataWeights ) throws Exception {
		initializeHiddenUniformly();
		//adjustParameters( classIndex, data, dataWeights, true, true, false );
		adjustParameters( classIndex, data, dataWeights, false, true, false );
		adjustParameters( classIndex, data, dataWeights, true, false, false );
	}
	
	protected void initializeUsingPlugIn( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception
	{
		double[] old = null;
		if( plugInBg ) {
			old = bg.getCurrentParameterValues();
		}
		bg.initializeFunction( index, freeParams, data, weights );

		// consensus 90%
		int num = getNumberOfMotifs(), a = (int) alphabets.getAlphabetLengthAt(0), s, p, l;
		double d = 0.1/(a-1), h;
		d = (1d-a*d)/(a*d);
		Sequence seq;
		for( int motif = 0; motif < num; motif++ ) {
			l = function[2*motif].getLength();
			s = r.nextInt( data[index].getNumberOfElements() );
			seq = data[index].getElementAt( s );
			p = r.nextInt( seq.getLength() - l + 1 );
			seq = seq.getSubSequence( p, l );
			h = d * function[2*motif].getESS();
			function[2*motif].initializeFunction( 0, freeParams, new DataSet[]{ new DataSet( "", seq ) }, new double[][]{{h}} );
			
			//System.out.print( s + "\t" + p + "\t" + seq + "\t" );
		}
		initializeHiddenUniformly();
		adjustParameters( index, data, weights, false, true, false ); // adjust duration
		adjustParameters( index, data, weights, true, false, true ); // adjust component parameters and motif
		if( plugInBg ) {
			bg.setParameters( old, 0 );
		} else {
			bg.initializeUniformly( freeParams );
		}
		//initBgHelp();
	}
	
	/**
	 * This method allows to adjust all parameter except those of the flanking sequence {@link DifferentiableStatisticalModel}.
	 * 
	 * @param data the {@link DataSet} containing sequences with motifs
	 * @param dataWeights the weights corresponding to the {@link Sequence}s in the {@link DataSet} <code>data</code>
	 * @param adjustComponent a switch to determine whether the component weights should be adjusted
	 * @param adjustDuration a switch to determine whether the {@link DurationDiffSM}s should be adjusted
	 * @param adjustMotif a switch to determine whether the {@link DifferentiableStatisticalModel}s that are used for the motifs should be adjusted
	 * 
	 * @throws Exception if some initialize or adjust method from {@link DifferentiableStatisticalModel} and {@link DurationDiffSM} is thrown 
	 */
	private void adjustParameters( int index, DataSet[] data, double[][] dataWeights, boolean adjustComponent, boolean adjustDuration, boolean adjustMotif ) throws Exception 
	{
		int c = getNumberOfComponents(), stop = c - (type == CONTAINS_ALWAYS_A_MOTIF?0:1), i = 0, n = 0, anz = data[index].getNumberOfElements(), l;

		int[][] len = null;
		DurationDiffSM dur;
		
		Sequence[][] seqs = null;
		double[][] seqComponentWeights = null;		
				
		
		double[][] weightsPos = null, weightsNeg = null;
		
		if( adjustDuration || adjustMotif ) {
			len = new int[stop][];
			for( i = 0; i < stop; i++ )
			{
				dur = (DurationDiffSM) function[2*i+1];
				len[i] = new int[dur.getNumberOfPossibilities()];
				dur.reset();
				l = 0;
				do
				{
					len[i][l++] = dur.internal[0];
				}
				while( dur.next() );
			}
			

			if( adjustMotif ) {
				seqs = new Sequence[getNumberOfMotifs()][anz];
				seqComponentWeights = new double[seqs.length][anz];
			}
		
			if( adjustDuration ) {
				weightsPos = new double[stop][];
				weightsNeg = new double[stop][];
				for( i = 0; i < stop; i++ )
				{
					weightsPos[i] = new double[len[i].length];
					weightsNeg[i] = new double[len[i].length];
				}
			}
		}
		double[] stat = new double[hiddenParameter.length];
		int maxIndex;
		double w = 1, scw;
		Sequence current;
		for( int d = 0; d < data.length; d++ )
		//int d = index;
		{
			anz = data[d].getNumberOfElements();
			for( n = 0; n < anz; n++ )
			{
				if( dataWeights != null )
				{
					w = dataWeights[d][n];
				}
				current = data[d].getElementAt( n );
				fillComponentScores( current, 0 );
				Normalisation.logSumNormalisation( componentScore );
				//if( adjustComponent ) System.out.println( Arrays.toString( componentScore ) );
				for( i = 0; i < c; i++ )
				{
					scw = componentScore[i]*w;

					if( adjustComponent && d == index ) {
						stat[i] += scw;
					}
					if( ( adjustDuration || adjustMotif ) && i < stop )
					{
						maxIndex = 0;
						for( l = 0; l < simpleScore[i].length; l++ )
						{
							if( simpleScore[i][l] >= simpleScore[i][maxIndex] )
							{
								maxIndex=l;
							}
						}
						
						if( adjustDuration ) {
							if( d == index ) {
								weightsPos[i][maxIndex] += scw;
							} else {
								weightsNeg[i][maxIndex] += scw;
							}							
						}
						
						if( adjustMotif && d == index ) {
							seqComponentWeights[i][n] = scw;
							seqs[i][n] = current.getSubSequence( len[i][maxIndex], function[2*i].getLength() );
							if( getStrandProbabilitiesFor( i, 0, seqs[i][n], 0 )[0] < 0.5 ) {
								seqs[i][n] = seqs[i][n].reverseComplement();
							}
						}
					}
				}
			}
		}
		
		//adjust parameters
		if( adjustComponent ) {
			computeHiddenParameter( stat, true );
		}
		if( adjustDuration ) {
			
			double f; //scaling factor to get pos. and neg. data on one scale
			double[] help;
			anz = data[index].getNumberOfElements();
			n = length / anz; //expected distance between two sites 
			c = Math.max( 1, n / 2 );
			double sumPos, sumNeg;
			/**/
			//TODO Problems with heuristic
			for( i = 0; i < stop; i++ ) {
				//System.out.println( Arrays.toString(weightsPos[i]) );
				//old: ((DurationDiffSM) function[2*i+1]).adjust( len[i], weightsPos[i] );
				/*
				String s = Arrays.toString( len[i] );
				System.out.println( "len=c(" + s.substring( 1, s.length()-1 ) + ");");
				s = Arrays.toString( weightsPos[i] );
				System.out.println( "wPos=c(" + s.substring( 1, s.length()-1 ) + ");");
				s = Arrays.toString( weightsNeg[i] );
				System.out.println( "wNeg=c(" + s.substring( 1, s.length()-1 ) + ");");/**/
				
				//smoothing of pos. and neg. weights (bins)
				help = new double[weightsPos[i].length];
				sumPos = 0;
				for( int binEnd, x, z, y = 0; y < help.length; y++ ) {
					help[y] = 0;
					binEnd = Math.min( help.length-1, y+c );
					for( z = 0, x = Math.max(0,y-c); x <= binEnd; x++, z++ ) {
						help[y] += weightsPos[i][x];
					}					
					help[y] /= z;
					sumPos += help[y];
				}
				weightsPos[i] = help;
				
				sumNeg = 0;
				for( int binEnd, x, z, y = 0; y < help.length; y++ ) {
					help[y] = 0;
					binEnd = Math.min( help.length, y+c );
					for( z = 0, x = Math.max(0,y-c); x < binEnd; x++, z++ ) {
						help[y] += weightsNeg[i][x];
					}					
					help[y] /= z;
					sumNeg += help[y];
				}
				weightsNeg[i] = help;
				/*
				s = Arrays.toString( weightsPos[i] );
				System.out.println( "wPos2=c(" + s.substring( 1, s.length()-1 ) + ");");
				s = Arrays.toString( weightsNeg[i] );
				System.out.println( "wNeg2=c(" + s.substring( 1, s.length()-1 ) + ");");/**/
				
				//contrast between pos and neg weights
				//System.out.println( sumPos + "\t" + sumNeg );
				f = sumNeg> 0 ? ( sumPos / sumNeg ) : 0;
				for( int x = 0; x < weightsPos.length; x++ ) {
					if( weightsPos[i][x] > 1.2*f*weightsNeg[i][x] ) {
						weightsPos[i][x] = weightsPos[i][x] - f*weightsNeg[i][x];
					} else {
						weightsPos[i][x] = 0;
					}
				}
				/*
				s = Arrays.toString( weightsPos[i] );
				System.out.println( "w=c(" + s.substring( 1, s.length()-1 ) + ");");
				System.out.println();

				/**/
				
				//System.out.println( Arrays.toString(weightsPos[i]) );
				((DurationDiffSM) function[2*i+1]).adjust( len[i], weightsPos[i] );
				//System.out.println( function[2*i+1] );
			}
		}
		if( adjustMotif ) {
			for( i = 0; i < stop; i++ ) {
				function[2*i].initializeFunction( 0, this.freeParams, new DataSet[]{new DataSet( "picked sites " + i, seqs[i] ) }, new double[][]{seqComponentWeights[i]} );
				//System.out.println( function[2*i+1] );
			}
		}
	}
	
	public double[] getStrandProbabilitiesFor( int component, int motif, Sequence sequence, int startpos ) throws Exception
	{	 
		int c = getNumberOfComponents() - (type==CONTAINS_SOMETIMES_A_MOTIF?1:0);
		if( motif > 0 || component >= c )
		{
			throw new IndexOutOfBoundsException();
		}
		else
		{
			DifferentiableSequenceScore m = function[2*component];
			while( m instanceof NormalizedDiffSM )
			{
				m = ((NormalizedDiffSM)m).getFunction();
			}
			if( m instanceof StrandDiffSM )
			{
				if(startpos == 0){
					return ((StrandDiffSM)m).getProbsForComponent( sequence );
				}else{
					return ((StrandDiffSM)m).getProbsForComponent( sequence.getSubSequence( startpos ) );
				}
			}
			else
			{
			    return new double[]{1.0,0.0};
			}
		}
	}
}
