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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

import de.jstacs.algorithms.graphs.UnionFind;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.DataSet.ElementEnumerator;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.CombinationIterator;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.InhCondProb;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.MEM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.MEMConstraint;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * The class manipulate some constraints.
 * 
 * @author Jens Keilwagen
 */
public class ConstraintManager
{
	private ConstraintManager()
	{
	}

	/**
	 * This method computes the (smoothed) relative frequencies. If <code>ess=0</code> no smoothing is done.
	 * @param ess
	 *            the ESS, if ESS is zero than MLE otherwise MAPE
	 * @param constr
	 *            the constraints, should be fill with absolute frequencies
	 * 
	 * @throws IllegalArgumentException
	 *             if the ess is negative
	 */
	public static void computeFreqs( double ess, Constraint... constr ) throws IllegalArgumentException
	{
		// check some constraints
		if( ess < 0 )
		{
			throw new IllegalArgumentException( "The ess has to be non-negative." );
		}

		// init the constraints with the pseudocounts
		for( int counter1 = 0; counter1 < constr.length; counter1++ )
		{
			constr[counter1].estimate( ess );
		}
	}

	/**
	 * Fills the (inhomogeneous) <code>constr</code> with the weighted absolute frequency of the DataSet
	 * <code>data</code> and computes the frequencies will not be computed.
	 * 
	 * @param alphabets
	 *            the alphabets over which the constraints are defined
	 * @param length
	 *            the length for which the constraints are defined
	 * @param data
	 *            the sequences
	 * @param weights
	 *            the weights for the sequences,
	 *            <ol>
	 *            <li> <code>weights==null</code> or
	 *            <li> <code>weights.length = data.getNumberOfElements()</code>, for all i:<code>weights[i]>=0</code>
	 *            </ol>
	 * @param reset whether the constraints should be reseted
	 * @param constr
	 *            constraints to fill
	 * @return the sum of the weights
	 * 
	 * @throws WrongAlphabetException if the alphabet of the data is not correct
	 * @throws IllegalArgumentException if the weights array has wrong dimension or the element length of the data is not correct
	 */
	public static double countInhomogeneous( AlphabetContainer alphabets, int length, DataSet data, double[] weights,
			boolean reset, Constraint... constr ) throws WrongAlphabetException, IllegalArgumentException
	{
		int d = data.getNumberOfElements(), counter1, counter2;
		Sequence seq;
		double all = 0;

		// check some constraints
		if( weights != null && d != weights.length )
		{
			throw new IllegalArgumentException( "The weights are not suitable for the data (wrong dimension)." );
		}
		if( !alphabets.checkConsistency( data.getAlphabetContainer() ) )
		{
			throw new WrongAlphabetException( "The alphabets of the model and the DataSet are not suitable." );
		}
		if( data.getElementLength() != length )
		{
			throw new IllegalArgumentException( "The sequence length of the model and the DataSet are not suitable." );
		}
		
		// reset counter
		if( reset )
		{
			for( counter1 = 0; counter1 < constr.length; counter1++ )
			{
				constr[counter1].reset();
			}
		}

		// fill the constraints with the absolute frequency in the data
		ElementEnumerator ei = new ElementEnumerator( data );
		if( weights == null )
		{
			all = d;
			for( counter1 = 0; counter1 < d; counter1++ )
			{
				seq = ei.nextElement();
				for( counter2 = 0; counter2 < constr.length; counter2++ )
				{
					constr[counter2].add( seq, 0, 1 );
				}
			}
		}
		else
		{
			for( counter1 = 0; counter1 < d; counter1++ )
			{
				seq = ei.nextElement();
				if( weights[counter1] > 0 )
				{
					all += weights[counter1];
					for( counter2 = 0; counter2 < constr.length; counter2++ )
					{
						constr[counter2].add( seq, 0, weights[counter1] );
					}
				}
				else if( weights[counter1] < 0 )
				{
					throw new IllegalArgumentException( "All weights have to be non-negative. Violated in position "
							+ counter1 + "." );
				}
			}
		}
		return all;
	}

	/**
	 * This method draws relative frequencies.
	 * @param ess
	 *            the ESS (additional pseudocount)
	 * @param constr
	 *            the constraints (can be fill with absolute frequencies)
	 * 
	 * @throws IllegalArgumentException
	 *             if the ess is negative
	 */
	public static void drawFreqs( double ess, InhCondProb... constr ) throws IllegalArgumentException
	{
		// check some constraints
		if( ess < 0 )
		{
			throw new IllegalArgumentException( "The ess has to be non-negative." );
		}

		// init the constraints with the pseudocounts
		for( int counter1 = 0; counter1 < constr.length; counter1++ )
		{
			constr[counter1].drawParameters( ess );
		}
	}
	
	/**
	 * Extracts the constraint of a String and returns an ArrayList of int[]. These can be used to create constraints.
	 * 
	 * @param length 
	 * 			  the sequence respectively the model length
	 * @param encoded
	 *            constraints encoded in a String
	 *            <ul>
	 *            	<li> items are separated by &quot;;&quot;
	 *            	<li> short notation for sets of constraints, e.g. &quot;m2sx&quot;
	 *            	<li> or each constraint as e.g. &quot;i,j,k;&quot; 
	 *            <ul>
	 * 
	 * @return an ArrayList of int[]
	 * 
	 * @throws IllegalArgumentException if the constraints can not be parsed correctly
	 */
	public static ArrayList<int[]> extract( int length, String encoded ) throws IllegalArgumentException
	{
		StringTokenizer simple, st = new StringTokenizer( encoded, ";" );
		String current, help;
		int mp, sp, m, i;
		ArrayList<int[]> list = new ArrayList<int[]>();
		int[] pos;
		boolean[][] constrAbbrev = new boolean[length + 1][];
		constrAbbrev[0] = null;
		constrAbbrev[1] = new boolean[]{ false };
		for( i = 2; i < constrAbbrev.length; i++ )
		{
			constrAbbrev[i] = new boolean[length - i + 1];
			Arrays.fill( constrAbbrev[i], false );
		}
		while( st.hasMoreTokens() )
		{
			current = st.nextToken().trim();
			mp = current.indexOf( "m" );
			sp = current.indexOf( "s" );
			if( mp >= 0 && sp > 0 )
			{
				// abbreviation for a set of constraints
				m = Integer.parseInt( current.substring( mp + 1, sp ) );
				if( m > length )
				{
					throw new IllegalArgumentException( "The marginal order of " + m + " is not possible for length "+ length );
				}
				help = current.substring( sp + 1 );
				if( help.equals( "x" ) )
				{
					for( i = 0; i < constrAbbrev[m].length; i++ )
					{
						constrAbbrev[m][i] = true;
					}
				}
				else
				{
					constrAbbrev[m][Integer.parseInt( current.substring( sp + 1 ) )] = true;
				}
			}
			else
			{
				// simple constraint
				m = 0;
				simple = new StringTokenizer( current, "," );
				pos = new int[simple.countTokens()];
				mp = 0;
				while( simple.hasMoreTokens() )
				{
					i = Integer.parseInt( simple.nextToken() );
					if( 0 <= i && i < length )
					{
						pos[mp++] = i;
					}
					else
					{
						throw new IllegalArgumentException( "Could not correctly parse \"" + current + "\"." );
					}
				}
				list.add( getSortedUnique( pos ) );
			}
		}
		add( constrAbbrev, list );
		return list;
	}

	private static String unfold( String constraint, int startpos, int endpos )
	{
		ArrayList<int[]> list = extract( endpos-startpos, constraint );
		StringBuffer erg = new StringBuffer( list.size()*10 );
		int[] current;
		for( int j, i = 0; i < list.size(); i++ )
		{
			current = list.get(i);
			erg.append( current[0]+startpos );
			for( j = 1; j < current.length; j++ )
			{
				erg.append( ","+(current[j]+startpos) );
			}
			erg.append( ";" );
		}
		return erg.toString();
	}

	/**
	 * Tries to compute the entropy as exact as possible.
	 * 
	 * @param c
	 *            the constraint
	 * 
	 * @return the entropy of the constraint
	 */
	public static double getEntropy( Constraint c )
	{
		double h = 0;
		int i = 0;
		double[] temps = new double[c.getNumberOfSpecificConstraints()];
		// 1. Schleife Summanden berechnen
		for( ; i < temps.length; i++ )
		{
			if( c.getFreq( i ) > 0 )
			{
				temps[i] = c.getFreq( i ) * Math.log( c.getFreq( i ) );
			}
		}

		// sortieren
		// => somit gleiche Reihenfolge (temps[0] <= temps[1] <= temps[2] <= ...)
		Arrays.sort( temps );

		// 2.Schleife Entropy aufaddieren
		for( i = 0; i < temps.length; i++ )
		{
			h -= temps[i];
		}
		return h;
	}

	/**
	 * Computes the sum of differences of the logarithmic values of the prior knowledge and all counts. <br>
	 * <br>
	 * In Latex notation:<br>
	 * <code>\sum_i [\ln \Gamma(\alpha_i) - \ln \Gamma(\alpha_i + N_i)]</code>
	 * 
	 * @param c
	 *            the constraint
	 * @param ess
	 *            the ESS
	 * 
	 * @return the sum
	 */
	public static double getLogGammaSum( Constraint c, double ess )
	{
		int i = 0, l = c.getNumberOfSpecificConstraints();
		double pc = ess / (double) l, sum = l * Gamma.logOfGamma( pc );
		while( i < l )
		{
			sum -= Gamma.logOfGamma( c.getCount( i++ ) + pc );
		}
		return sum;
	}

	/**
	 * This method tries to find and remove subconstraints that are already fulfilled by a bigger one. The method only
	 * looks for the positions, so it is recommended to use this method before any learning step.
	 * 
	 * @param list
	 *            the list of all constraints
	 */
	public static void reduce( AbstractList<int[]> list )
	{
		int counter1 = 0, counter2;
		// find subconstraints
		int[] current, helpPos;
		while( counter1 < list.size() )
		{
			helpPos = list.get( counter1 );
			counter2 = counter1 + 1;
			while( counter2 < list.size() )
			{
				current = list.get( counter2 );
				if( helpPos.length < current.length )
				{
					if( isSubset( helpPos, current ) )
					{
						list.remove( counter1 );
						counter1--;
						break;
					}
					else
					{
						counter2++;
					}
				}
				else
				{
					if( isSubset( current, helpPos ) )
					{
						list.remove( counter2 );
					}
					else
					{
						counter2++;
					}
				}
			}
			counter1++;
		}
	}

	private static void add( boolean[][] constrAbbrev, ArrayList<int[]> list )
	{
		int length = constrAbbrev.length - 1, max = 0;
		int counter1 = 1, counter2=1, distance;
		// find maximal MEM-order that will be used 
		for( ; counter2 <= length; counter2++ )
		{
			for( distance = 0; distance < constrAbbrev[counter2].length; distance++ )
			{
				if( constrAbbrev[counter2][distance] )
				{
					max = counter2;
				}
			}
		}
		CombinationIterator c = new CombinationIterator( length, max );
		int[] current;
		for( ; counter1 <= max; counter1++ )
		{
			counter2 = 0;
			while( counter2 < constrAbbrev[counter1].length && !constrAbbrev[counter1][counter2] )
			{
				counter2++;
			}
			if( counter2 != constrAbbrev[counter1].length )
			{
				c.setCurrentLength( counter1 );
				do
				{
					current = c.getCombination();
					distance = 1 - counter1;
					counter2 = 1;
					while( counter2 < counter1 )
					{
						distance += current[counter2] - current[counter2++ - 1];
					}
					if( constrAbbrev[counter1][distance] )
					{
						list.add( current );
					}
				}
				while( c.next() );
			}
		}
	}

	/**
	 * Returns the intersection of the sorted arrays <code>array1</code> and <code>array2</code>.
	 */
	private static int[] getIntersection( int[] array1, int[] array2 )
	{
		int counter1 = 0, counter2 = 0, length = 0;
		int[] erg, zerg = new int[Math.min( array1.length, array2.length )];
		while( counter1 < array1.length && counter2 < array2.length )
		{
			// find next
			while( counter1 < array1.length && array1[counter1] < array2[counter2] )
			{
				counter1++;
			}
			if( counter1 < array1.length )
			{
				while( counter2 < array2.length && array1[counter1] > array2[counter2] )
				{
					counter2++;
				}
				// set common value as part of the result
				if( counter2 < array2.length && array1[counter1] == array2[counter2] )
				{
					zerg[length++] = array1[counter1++];
					counter2++;
				}
			}
		}
		erg = new int[length];
		System.arraycopy( zerg, 0, erg, 0, length );
		return erg;
	}

	/**
	 * Returns a sorted array with unique entries.
	 * 
	 * @param pos
	 *            the array
	 * 
	 * @return the sorted array
	 * 
	 * @throws IllegalArgumentException
	 *             if the values of pso are not unique
	 */
	private static int[] getSortedUnique( int[] pos ) throws IllegalArgumentException
	{
		int[] sorted = new int[pos.length];
		System.arraycopy( pos, 0, sorted, 0, pos.length );
		Arrays.sort( sorted );
		int i = 1;
		while( i < pos.length && sorted[i - 1] < sorted[i] )
		{
			i++;
		}
		if( i < sorted.length )
		{
			throw new IllegalArgumentException( "The position array is not unique." );
		}
		return sorted;
	}

	/**
	 * <b>The &quot;sets&quot; have to be sorted and unique.</b>
	 */
	private static boolean isSubset( int[] candidateSubset, int[] set )
	{
		int counter1 = 0, counter2 = 0;
		while( counter2 < candidateSubset.length )
		{
			while( counter1 < set.length && set[counter1] < candidateSubset[counter2] )
			{
				counter1++;
			}
			if( counter1 == set.length || set[counter1] > candidateSubset[counter2] )
			{
				return false;
			}
			else
			{
				counter2++;
				counter1++;
			}
		}
		return true;
	}

	private static ArrayList[] split( AbstractList<int[]> list, int length, int[][] indices )
	{
		int counter1 = 0, counter2;
		ArrayList[] constr = new ArrayList[indices.length];
		int[] current, component = new int[length];
		for( ; counter1 < indices.length; counter1++ )
		{
			constr[counter1] = new ArrayList<int[]>();
			for( counter2 = 0; counter2 < indices[counter1].length; counter2++ )
			{
				component[indices[counter1][counter2]] = counter1;
			}
		}
		for( counter1 = 0; counter1 < list.size(); counter1++ )
		{
			current = list.get( counter1 );
			constr[component[current[0]]].add( current );
		}
		return constr;
	}

	private static int[][] unionFind( AbstractList<int[]> list, int[] alphabetLength, int leafOutIndex )
	{
		UnionFind uf = new UnionFind( alphabetLength.length );
		int[] current;
		int counter1 = 0, counter2;
		for( ; counter1 < list.size(); counter1++ )
		{
			if( counter1 != leafOutIndex )
			{
				current = list.get( counter1 );
				for( counter2 = 1; counter2 < current.length; counter2++ )
				{
					uf.union( current[0], current[counter2] );
				}
			}
		}
		return uf.getComponents();
	}
	
	/**
	 * Creates the constraints for a part of a model
	 * 
	 * @param list
	 *            the list of all cliques, each clique is used for one constraint
	 * @param alphabetLength
	 *            the array of alpahebtLength for each position
	 * @param indices
	 *            the positions used in this part of the model
	 * 
	 * @return the array of constraints
	 */
	public static MEMConstraint[] createConstraints( AbstractList<int[]> list, int[] alphabetLength, int[] indices )
	{
		if( alphabetLength.length == indices.length )
		{
			return createConstraints( list, alphabetLength );
		}
		else
		{
			int[] corrected = new int[alphabetLength.length];
			int i = 0, j;
			while( i < indices.length )
			{
				corrected[indices[i]] = i++;
			}
			int[] pos, corPos;
			MEMConstraint[] constr = new MEMConstraint[list.size()];
			for( i = 0; i < constr.length; i++ )
			{
				pos = list.get( i );
				corPos = new int[pos.length];
				for( j = 0; j < pos.length; j++ )
				{
					corPos[j] = corrected[pos[j]];
				}
				constr[i] = new MEMConstraint( pos, alphabetLength, corPos );
			}
			return constr;
		}
	}
	
	/**
	 * Creates the constraints of a model
	 * 
	 * @param list
	 *            the list of all cliques, each clique is used for one constraint
	 * @param alphabetLength
	 *            the array of alpahebtLength for each position
	 * 
	 * @return the array of constraints
	 */
	public static MEMConstraint[] createConstraints( AbstractList<int[]> list, int[] alphabetLength )
	{
		MEMConstraint[] constr = new MEMConstraint[list.size()];
		for( int i = 0; i < constr.length; i++ )
		{
			constr[i] = new MEMConstraint( list.get( i ), alphabetLength );
		}
		return constr;
	}
	
	private static void addDivided( ArrayList<MEM> mems, ArrayList<int[]> constr, int[] indices, int[] alphabetLength )
	{
		if( constr.size() == 1 )
		{
			mems.add( new MEM( constr.get( 0 ), alphabetLength, null ) );
		}
		else
		{
			int[][] parts;
			int[] current;
			ArrayList<int[]> cond;
			int counter1 = constr.size() - 1, counter2;
			while( counter1 >= 0 )
			{
				parts = unionFind( constr, alphabetLength, counter1 );
				if( parts.length > 1 )
				{
					current = constr.remove( counter1 );
					cond = new ArrayList<int[]>();
					ArrayList[] splittedConstr = split( constr, alphabetLength.length, parts );
					// parts.length == splittedConstr.length()
					for( counter2 = 0; counter2 < parts.length; counter2++ )
					{
						if( splittedConstr[counter2].size() > 0 )
						{
							cond.add( getIntersection( current, parts[counter2] ) );
						}
					}
					// add a new MEM for the edge left out
					mems.add( new MEM( current, alphabetLength, cond.toArray( new int[0][] ) ) );
					for( counter2 = 0; counter2 < parts.length; counter2++ )
					{
						if( splittedConstr[counter2].size() > 0 )
						{
							addDivided( mems, (ArrayList<int[]>) splittedConstr[counter2], parts[counter2],
									alphabetLength );
						}
					}
					// all work is done
					return;
				}
				else
				{
					// nothing to do, step to next constraint
					counter1--;
				}
			}
			mems.add( new MEM( constr, alphabetLength, indices ) );
		}
	}
	
	
	/**
	 * This enum defines the different possible types of decomposition of a model.
	 * 
	 * @author Jens Keilwagen
	 */
	public static enum Decomposition {
		/**
		 * This value indicates that the model should have only one component.
		 */
		DECOMPOSE_NOTHING,
	
		/**
		 * This value indicates that the model should be decomposed in its unconnected components.
		 */
		DECOMPOSE_UNCONNECTED,
	
		/**
		 * This value indicates that the model should be decomposed in components that were connected, but
		 * could be separated.
		 */
		DECOMPOSE_LESS_CONNECTED;
	}
	
	/**
	 * This method tries to disconnect the constraints and create the models. This can be useful to be much faster
	 * while training.
	 * 
	 * @param list
	 *            the list of all constraints
	 * @param alphabetLength
	 *            the length of the alphabets at the specific positions
	 * @param decomposition
	 *            the choice how to disconnect the constraints
	 *            
	 * @return an array of {@link MEM}
	 * 
	 * @see Decomposition#DECOMPOSE_LESS_CONNECTED
	 * @see Decomposition#DECOMPOSE_UNCONNECTED
	 * @see Decomposition#DECOMPOSE_NOTHING
	 */
	public static MEM[] disconnect( AbstractList<int[]> list, int[] alphabetLength, Decomposition decomposition )
	{
		if( decomposition == Decomposition.DECOMPOSE_NOTHING )
		{
			// all positions are in one component
			int[] ind = new int[alphabetLength.length];
			for( int counter1 = 0; counter1 < ind.length; counter1++ )
			{
				ind[counter1] = counter1;
			}
			// make one component
			return new MEM[]{ new MEM( list, alphabetLength, ind ) };
		}
		else if( decomposition == Decomposition.DECOMPOSE_UNCONNECTED || decomposition == Decomposition.DECOMPOSE_LESS_CONNECTED )
		{
			// find unconnected components by using union find
			int[][] indices = unionFind( list, alphabetLength, -1 );
			// make as many components as you can find with union find
			ArrayList[] constr = split( list, alphabetLength.length, indices );
			if( decomposition == Decomposition.DECOMPOSE_UNCONNECTED )
			{
				MEM[] res = new MEM[indices.length];
				for( int counter1 = 0; counter1 < indices.length; counter1++ )
				{
					res[counter1] = new MEM( (ArrayList<int[]>) constr[counter1], alphabetLength, indices[counter1] );
				}
				return res;
			}
			else
			{
				ArrayList<MEM> mems = new ArrayList<MEM>( indices.length * 2 );
				for( int counter1 = 0; counter1 < indices.length; counter1++ )
				{
					addDivided( mems, (ArrayList<int[]>) constr[counter1], indices[counter1], alphabetLength );
				}
				return mems.toArray( new MEM[0] );
			}
		}
		else
		{
			throw new IllegalArgumentException( "The choice of decomposition was illegal." );
		}
	}
}
